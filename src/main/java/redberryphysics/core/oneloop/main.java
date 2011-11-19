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
package redberryphysics.core.oneloop;

import java.util.ArrayList;
import java.util.List;
import redberry.core.context.CC;
import redberry.core.context.ToStringMode;
import redberry.core.tensor.Expression;
import redberry.core.tensor.MultiTensor;
import redberry.core.tensor.Product;
import redberry.core.tensor.SimpleTensor;
import redberry.core.tensor.Sum;
import redberry.core.tensor.Tensor;
import redberry.core.tensor.TensorIterator;
import redberry.core.tensor.TensorNumber;
import redberry.core.tensor.iterators.TensorFirstTreeIterator;
import redberry.core.tensor.iterators.TensorTreeIterator;
import redberry.core.tensor.test.TTest;
import redberry.core.transformation.CalculateNumbers;
import redberry.core.transformation.Transformation;
import redberry.core.transformation.collect.EqualsSplitCriteria;
import redberry.core.transformation.collect.ScalarsSplitCriteria;
import redberry.core.transformation.concurrent.ExpandAndCollectTransformation;
import redberry.core.transformation.contractions.IndexesContractionsTransformation;
import redberry.core.utils.Indicator;
import redberryphysics.core.util.SqrSubs;

/**
 *
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
public class main {
    public static void main(String[] args) {
        OneLoop loop1 = new OneLoop(OneLoop.EVAL.INITIALIZE);
        loop1.evalRR();
        loop1.evalHatK();
        loop1.evalDeltas();
        for (Expression delta : loop1.DELTAs)
            loop1.RR.eval(delta.asSubstitution());
        for (Expression hatK : loop1.HATKs)
            loop1.RR.eval(hatK.asSubstitution());
        loop1.RR.eval(
                IndexesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                loop1.KRONECKER_DIMENSION.asSubstitution());

        Tensor rhs = ((MultiTensor) loop1.RR.right()).getElements().get(0).clone();
        System.out.println(" Evaluating RR ");

        TensorTreeIterator treeIterator = new TensorFirstTreeIterator(rhs);
        Tensor current;

        List<Sum> sums = new ArrayList<>();
        TensorIterator iterator = rhs.iterator();
        while (iterator.hasNext()) {
            current = iterator.next();
            if (current instanceof Sum && !TTest.testIsSymbol(current)) {
                sums.add((Sum) current.clone());
                iterator.remove();
            }
        }

        System.out.println(" counting terms: ... Sums number " + sums.size());
        Transformation tr = new SqrSubs((SimpleTensor) CC.parse("n_{\\alpha}"));

        long start = System.currentTimeMillis();
        Transformation ec = new ExpandAndCollectTransformation(
                EqualsSplitCriteria.INSTANCE1,
                Indicator.SYMBOL_INDICATOR,
                new Transformation[]{
                    IndexesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                    loop1.KRONECKER_DIMENSION.asSubstitution(),
                    tr,
                    CalculateNumbers.INSTANCE});


        Transformation tr1 = new Transformation() {
            @Override
            public Tensor transform(Tensor tensor) {
                if (!(tensor instanceof Product))
                    return tensor;
                TensorIterator iterator = tensor.iterator();
                while (iterator.hasNext())
                    if (TTest.testIsSymbol(iterator.next()))
                        iterator.remove();
                return tensor.equivalent();
            }
        };

        iterator = sums.get(1).iterator();
        while (iterator.hasNext())
            iterator.set(tr1.transform(iterator.next()));
        iterator = sums.get(2).iterator();
        while (iterator.hasNext())
            iterator.set(tr1.transform(iterator.next()));

        Product firstPart = new Product(sums.get(1), sums.get(2));
        System.out.println(sums.get(1).size());
        System.out.println(sums.get(2).size());
        System.out.println(sums.get(1).toString(ToStringMode.REDBERRY_SOUT));
        System.out.println(sums.get(2).toString(ToStringMode.REDBERRY_SOUT));
        Tensor fst = ec.transform(firstPart);
        System.out.println("+");
//        Tensor scn = ec.transform(secondPart);
        System.out.println(((MultiTensor) fst).size());
//        System.out.println(((MultiTensor) scn).size());
        Product fin = new Product();
        fin.add(rhs);
        fin.add(fst);
        fin.add(sums.get(0));
//        
        Transformation sP = new Transformation() {
            @Override
            public Tensor transform(Tensor tensor) {
                if (!(tensor instanceof Product))
                    throw new RuntimeException();
                Product s = new Product();
                for (Tensor t : tensor)
                    if (TTest.testIsSymbol(t))
                        s.add(t);
                if (s.isEmpty())
                    return TensorNumber.createONE();
                else
                    return s.equivalent();
            }
        };
        Transformation ec1 = new ExpandAndCollectTransformation(
                ScalarsSplitCriteria.INSTANCE,
                Indicator.SYMBOL_INDICATOR,
                new Transformation[]{
                    IndexesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                    loop1.KRONECKER_DIMENSION.asSubstitution(),
                    CalculateNumbers.INSTANCE,
                    sP});

        Tensor AAA = ec1.transform(rhs);
        System.out.println(((MultiTensor) AAA).size());
//        System.out.println(AAA.toString(ToStringMode.UTF8));
//        for (Sum sum : sums) {
//            rhs = new Product(rhs.equivalent(), sum);
//            rhs = ec.transform(rhs.clone());
//            System.out.println("+");
//        }

//        Transformation ec2 = new ExpandAndCollectTransformation(
//                ScalarsSplitCriteria.INSTANCE,
//                Indicator.FALSE_INDICATOR,
//                new Transformation[]{
//                    CalculateNumbers.INSTANCE});
//        Tensor Res = ec2.transform(AAA);
//        System.out.println(((MultiTensor) Res).size());
//        System.out.println(Res);
//        long stop = System.currentTimeMillis();
//        System.out.println("... . Evalution time = " + (stop - start) + "ms");
//        System.out.println(rhs.toString(ToStringMode.UTF8));

    }
}
