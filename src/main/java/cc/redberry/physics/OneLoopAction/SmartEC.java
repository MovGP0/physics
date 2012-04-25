/*
 * Redberry: symbolic tensor computations.
 *
 * Copyright (c) 2010-2012:
 *   Stanislav Poslavsky   <stvlpos@mail.ru>
 *   Bolotin Dmitriy       <bolotin.dmitriy@gmail.com>
 *
 * This file is part of Redberry.
 *
 * Redberry is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Redberry is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Redberry. If not, see <http://www.gnu.org/licenses/>.
 */
package cc.redberry.physics.OneLoopAction;

import cc.redberry.core.tensor.Product;
import cc.redberry.core.tensor.Sum;
import cc.redberry.core.tensor.Tensor;
import cc.redberry.core.tensor.TensorIterator;
import cc.redberry.core.tensor.testing.TTest;
import cc.redberry.transformation.Transformation;
import cc.redberry.transformation.Transformations;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
public class SmartEC implements Transformation {
    final Transformation ec;

    public SmartEC(Transformation ec) {
        this.ec = ec;
    }

    @Override
    public Tensor transform(Tensor tensor) {
          if (tensor instanceof Sum) {
            Sum sum = (Sum) tensor;
            List<Tensor> result = new ArrayList<>();
            int i = 0;
            for (Tensor summand : sum)
//                System.out.println(summand.toString(ToStringMode.UTF8));
                result.add(transform(summand)); //                System.out.println(result.get(i++));
            return new Sum(result);
        }
        if (!(tensor instanceof Product))
            return tensor;
        ArrayList<Tensor> sums = new ArrayList<>();
        TensorIterator iterator = tensor.iterator();
        while (iterator.hasNext()) {
            Tensor t = iterator.next().equivalent();
            if (t instanceof Sum && !TTest.testIsSymbol(t)) {
                sums.add(t);
                iterator.remove();
            }
        }
        ArrayList<Tensor> sumsNext = new ArrayList<>(), sumsTmp;
        while (sums.size() > 1) {
            for (int i = 0; i < sums.size() / 2; ++i)
                if (i * 2 + 1 < sums.size()) {
//                    System.out.print("Iter: " + ((Sum) sums.get(i * 2)).size() + ", " + ((Sum) sums.get(i * 2 + 1)).size() + "; ");
                    Tensor t = Transformations.expandAndCollectAllScalars(
                            ec.transform(new Product(sums.get(i * 2), sums.get(i * 2 + 1))));
                    //TODO review
                    if (t instanceof Sum)
                        sumsNext.add(
                                Transformations.expandAndCollectAllScalars(
                                ec.transform(new Product(sums.get(i * 2), sums.get(i * 2 + 1)))));
                    else
                        ((Product) tensor).add(t);
                }
            if (sums.size() % 2 == 1)
                sumsNext.add(sums.get(sums.size() - 1));
            sumsTmp = sums;
            sums = sumsNext;
            sumsNext = sumsTmp;
            sumsNext.clear();
//            System.out.println("Done.");
        }
        //TODO review
        if (sums.isEmpty())
            return ec.transform(tensor);
        if (((Product) tensor).isEmpty())
            return sums.get(0);
        ((Product) tensor).add(sums.get(0));
        return ec.transform(tensor);
    }
}
