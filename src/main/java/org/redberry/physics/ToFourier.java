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
package org.redberry.physics;

import java.util.ArrayList;
import java.util.List;
import redberry.core.context.CC;
import redberry.core.indexes.IndexesFactory;
import redberry.core.indexes.SimpleIndexes;
import redberry.core.indexgenerator.IntGenerator;
import redberry.core.tensor.Derivative;
import redberry.core.tensor.Integral;
import redberry.core.tensor.Product;
import redberry.core.tensor.SimpleTensor;
import redberry.core.tensor.Sum;
import redberry.core.tensor.Tensor;
import redberry.core.tensor.TensorField;
import redberry.core.tensor.TensorIterator;
import redberry.core.tensor.TensorNumber;
import redberry.core.tensor.indexmapping.IndexMappingDirect;
import redberry.core.tensor.iterators.TensorFirstTreeIterator;
import redberry.core.tensor.iterators.TensorTreeIterator;
import redberry.core.transformation.ApplyIndexMappingDirectTransformation;
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
        if (tensor instanceof Derivative) {
            Derivative d = (Derivative) tensor;
            Tensor[] vars = d.getVars();

            boolean c = true;
            for (int i = 0; i < vars.length; ++i)
                if (((SimpleTensor) vars[i]).getName() != from.getName()) {
                    c = false;
                    break;
                }
            if (c)
                if (mode(d.getTarget()) == TransformMode.TRANSFORM)
                    return TransformMode.TRANSFORM;
                else
                    return TransformMode.ERROR;
        }
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
                    throw new UnsupportedOperationException("Not supported for Dirac delta.");
            }
        }
        if (target instanceof Derivative)
            switch (mode(target)) {
                case TRANSFORM:
                    return TensorNumber.createZERO();
                case SKIP:
                case ERROR:
                    throw new UnsupportedOperationException("Not supported for Dirac delta.");
            }
        if (target instanceof Product) {
            final Product product = (Product) target.clone();
            final TensorIterator iterator = product.iterator();
            final List<Tensor> toFourier = new ArrayList<>();
            Tensor current;
            out_while:
            while (iterator.hasNext()) {
                current = iterator.next();
                switch (mode(current)) {
                    case ERROR:
                        throw new UnsupportedOperationException("Not supported.");
                    case SKIP:
                        continue out_while;
                    case TRANSFORM:
                        toFourier.add(current);
                        iterator.remove();
                }
            }

            if (toFourier.isEmpty())
                throw new UnsupportedOperationException("Not supported.");
            if (toFourier.size() == 1) {
                final Tensor tf = toFourier.get(0);
                if (tf instanceof Derivative)
                    return TensorNumber.createZERO();
                if (tf instanceof TensorField)
                    product.add(CC.createTensorField(((TensorField) tf).getName(), ((TensorField) tf).getIndexes().clone(), TensorNumber.createZERO()));
                return product;
            }

            final Generator generator = new Generator(product, to);
            final int size = toFourier.size();
            final List<SimpleTensor> generatedArgs = new ArrayList<>();
            SimpleTensor currentArg;
            for (int i = 0; i < size; ++i) {
                if (i == size - 1) {
                    final Sum sum = new Sum();
                    for (Tensor generated : generatedArgs)
                        sum.add(new Product(TensorNumber.createMINUSONE(), generated.clone()));
                    product.add(processTensor(toFourier.get(i), sum.equivalent(), generatedArgs.get(0).getIndexes()));
                    return new Integral(product, generatedArgs.toArray(new SimpleTensor[generatedArgs.size()]));
                }
                currentArg = generator.getNext();
                generatedArgs.add(currentArg);
                product.add(processTensor(toFourier.get(i), currentArg, currentArg.getIndexes()));
            }
        }
        throw new UnsupportedOperationException("Not supported.");
    }

    private Tensor processTensor(final Tensor t, final Tensor newArg, final SimpleIndexes newArgSimpleIndexes) {
        if (t instanceof TensorField) {
            final TensorField field = (TensorField) t;
            return CC.createTensorField(field.getName(), IndexesFactory.createSimple(field.getIndexes()), newArg.clone());
        }
        if (t instanceof Derivative) {
            final Derivative d = (Derivative) t;
            final Tensor[] vars = d.getVars();
            final int vC = vars.length;
            final TensorNumber image =
                    vC % 4 == 0 ? TensorNumber.createONE()
                    : (vC + 1) % 4 == 0 ? TensorNumber.createIMAGE_MINUSONE()
                    : vC % 2 == 0 ? TensorNumber.createMINUSONE() : TensorNumber.createIMAGE_ONE();
            final Product result = new Product();
            result.add(image);
            IndexMappingDirect im;
            Tensor newArg_;
            for (int i = 0; i < vC; ++i) {
                newArg_ = newArg.clone();
                im = new IndexMappingDirect(newArgSimpleIndexes, vars[i].getIndexes());
                newArg_ = ApplyIndexMappingDirectTransformation.INSTANCE.perform(newArg_, im);
                result.add(newArg_);
            }
            result.add(processTensor(d.getTarget(), newArg.clone(), newArgSimpleIndexes));
            return result;
        }
        throw new RuntimeException();
    }

    private class Generator {
        private final String sampleName;
        private final IntGenerator generator;
        private final SimpleIndexes sampleIndexes;

        Generator(final Tensor tensor, final SimpleTensor sample) {
            sampleIndexes = sample.getIndexes();
            sampleName = CC.getNameDescriptor((sample).getName()).getName();
            final List<Integer> integers = new ArrayList<>();
            final TensorTreeIterator iterator = new TensorFirstTreeIterator(tensor);
            Tensor current;
            while (iterator.hasNext()) {
                current = iterator.next();
                if (!(current instanceof SimpleTensor))
                    continue;
                final String name = CC.getNameDescriptor(((SimpleTensor) current).getName()).getName();

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
                final int[] engagedData = new int[integers.size()];
                for (int i = 0; i < integers.size(); ++i)
                    engagedData[i] = integers.get(i).intValue();
                generator = new IntGenerator(engagedData);
            }
        }
        private boolean first = true;

        SimpleTensor getNext() {
            if (first) {
                first = false;
                return CC.createSimpleTensor(sampleName, IndexesFactory.createSimple(sampleIndexes));
            } else
                return CC.createSimpleTensor(sampleName + generator.getNext(), IndexesFactory.createSimple(sampleIndexes));
        }
    }
}
