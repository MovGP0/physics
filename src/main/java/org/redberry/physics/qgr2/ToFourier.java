/*
 *
 * Redberry: symbolic tensor computations library.
 * Copyright (C) 2010-2011  Stanislav Poslavsky <stvlpos@mail.ru>
 *
 * This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
package org.redberry.physics.qgr2;

import java.util.ArrayList;
import java.util.List;
import redberry.core.context.CC;
import redberry.core.indexes.SimpleIndexes;
import redberry.core.indexgenerator.IntGenerator;
import redberry.core.tensor.Integral;
import redberry.core.tensor.Product;
import redberry.core.tensor.SimpleTensor;
import redberry.core.tensor.Sum;
import redberry.core.tensor.Tensor;
import redberry.core.tensor.TensorField;
import redberry.core.tensor.TensorIterator;
import redberry.core.tensor.TensorNumber;
import redberry.core.tensor.iterators.TensorFirstTreeIterator;
import redberry.core.tensor.iterators.TensorTreeIterator;
import redberry.core.transformation.Transformation;
import redberry.core.utils.TensorUtils;

/**
 *
 * @author Stanislav Poslavsky
 */
public class ToFourier implements Transformation {
    private final SimpleTensor from, to;

    public ToFourier(final SimpleTensor from, final SimpleTensor to) {
        this.from = from;
        this.to = to;
    }

    public ToFourier(final String from, final String to) {
        this((SimpleTensor) CC.parse(from), (SimpleTensor) CC.parse(to));
    }

    private boolean valid(final Integral integral) {
        final SimpleTensor[] vars = integral.vars();
        if (vars.length != 1)
            return false;
        if (vars[0].getName() != from.getName())
            return false;
        return true;
    }

    private static enum TransformMode {
        TRANSFORM, SKIP, ERROR;
    }

    private TransformMode mode(final Tensor tensor) {
        if (tensor instanceof SimpleTensor) {
            final int name = ((SimpleTensor) tensor).getName();
            if (name == from.getName() || name == to.getName())
                return TransformMode.ERROR;
            if (tensor.getClass() == SimpleTensor.class)
                return TransformMode.SKIP;
            if (!(tensor instanceof TensorField))
                return TransformMode.ERROR;

            final Tensor[] args = ((TensorField) tensor).getArgs();
            if (args.length == 1) {
                final Tensor arg = args[0];
                if (arg instanceof SimpleTensor)
                    if (((SimpleTensor) arg).getName() == from.getName())
                        return TransformMode.TRANSFORM;
                    else if (((SimpleTensor) arg).getName() == to.getName())
                        return TransformMode.ERROR;
                    else
                        return TransformMode.SKIP;
                if (TensorUtils.contains(arg, from, to))
                    return TransformMode.ERROR;
            }
            for (int i = 0; i < args.length; ++i)
                if (TensorUtils.contains(args[i], from, to))
                    return TransformMode.ERROR;
            return TransformMode.SKIP;
        }
        if (TensorUtils.contains(tensor, from, to))
            return TransformMode.ERROR;
        return TransformMode.SKIP;
    }

    @Override
    public Tensor transform(final Tensor tensor) {
        if (!(tensor instanceof Integral))
            return tensor;
        final Integral integral = (Integral) tensor;
        if (!valid(integral))
            throw new UnsupportedOperationException("Not supported.");
        final Tensor target = integral.target();
        if (target instanceof SimpleTensor) {
            final SimpleTensor st = (SimpleTensor) target;
            switch (mode(st)) {
                case TRANSFORM:
                    return CC.createTensorField(st.getName(), st.getIndexes(), TensorNumber.createZERO());
                case SKIP:
                case ERROR:
                    throw new UnsupportedOperationException("Not supported.");
            }
        }
        if (target instanceof Product) {
            final Product product = (Product) target.clone();
            final TensorIterator iterator = product.iterator();
            final List<TensorField> toFourier = new ArrayList<>();
            Tensor current;
            while (iterator.hasNext()) {
                current = iterator.next();
                switch (mode(current)) {
                    case ERROR:
                        throw new UnsupportedOperationException("Not supported.");
                    case SKIP:
                        continue;
                    case TRANSFORM:
                        toFourier.add((TensorField) current);
                        iterator.remove();
                }
            }

            if (toFourier.isEmpty())
                throw new UnsupportedOperationException("Not supported.");
            if (toFourier.size() == 1) {
                final TensorField tf = toFourier.get(0);
                product.add(CC.createTensorField(tf.getName(), tf.getIndexes(), TensorNumber.createZERO()));
                return product;
            }

            final Generator generator = new Generator(product, to);
            final int size = toFourier.size();
            final List<SimpleTensor> generatedArgs = new ArrayList<>();
            SimpleTensor currentArg;
            for (int i = 0; i < size; ++i) {
                if (i == size - 1) {
                    final Sum sum = new Sum();
                    sum.add(to);
                    for (Tensor generated : generatedArgs)
                        sum.add(new Product(TensorNumber.createMINUSONE(), generated));
                    product.add(CC.createTensorField(to.getName(), to.getIndexes(), sum.equivalent()));
                    return new Integral(product, generatedArgs.toArray(new SimpleTensor[generatedArgs.size()]));
                }
                currentArg = generator.getNext();
                generatedArgs.add(currentArg);
                product.add(CC.createTensorField(to.getName(), to.getIndexes(), currentArg));
            }
        }
        throw new UnsupportedOperationException("Not supported.");
    }

    private class Generator {
        final String sampleName;
        final IntGenerator generator;
        final SimpleIndexes indexes;

        Generator(final Tensor tensor, final SimpleTensor sample) {
            indexes = sample.getIndexes();
            sampleName = CC.getNameDescriptor((sample).getName()).getName();
            final List<Integer> integers = new ArrayList<>();
            final TensorTreeIterator iterator = new TensorFirstTreeIterator(tensor);
            Tensor current;
            while (iterator.hasNext()) {
                current = iterator.next();
                if (!(current instanceof SimpleTensor))
                    continue;
                String name = CC.getNameDescriptor(((SimpleTensor) current).getName()).getName();

                if (name.length() < sampleName.length() || !name.substring(0, sampleName.length()).equals(sampleName))
                    continue;
                Integer value;
                try {
                    value = Integer.parseInt(name.substring(sampleName.length()));
                } catch (NumberFormatException exception) {
                    continue;
                }
                integers.add(value);
            }
            if (integers.isEmpty())
                generator = new IntGenerator();
            else {
                int[] engagedData = new int[integers.size()];
                for (int i = 0; i < integers.size(); ++i)
                    engagedData[i] = integers.get(i).intValue();
                generator = new IntGenerator(engagedData);
            }
        }

        SimpleTensor getNext() {
            int val = generator.getNext();
            return CC.createSimpleTensor(sampleName + val, indexes);
        }
    }
}
