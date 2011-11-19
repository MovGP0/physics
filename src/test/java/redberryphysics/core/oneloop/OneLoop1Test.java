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

import redberry.core.tensor.SimpleTensor;
import redberryphysics.core.util.SqrSubs;
import java.util.ArrayList;
import java.util.List;
import org.junit.Ignore;
import redberry.core.tensor.Sum;
import redberry.core.tensor.Expression;
import redberry.core.transformation.Transformation;
import redberry.core.tensor.Tensor;
import org.junit.Test;
import redberry.core.context.CC;
import redberry.core.context.ToStringMode;
import redberry.core.tensor.MultiTensor;
import redberry.core.tensor.Product;
import redberry.core.tensor.TensorIterator;
import redberry.core.tensor.TensorNumber;
import redberry.core.tensor.iterators.TensorFirstTreeIterator;
import redberry.core.tensor.iterators.TensorTreeIterator;
import redberry.core.tensor.test.TTest;
import redberry.core.transformation.CalculateNumbers;
import redberry.core.transformation.ExpandBrackets;
import redberry.core.transformation.IndexesInsertion;
import redberry.core.transformation.Transformations;
import redberry.core.transformation.Transformer;
import redberry.core.transformation.collect.CollectFactory;
import redberry.core.transformation.collect.CollectInputPort;
import redberry.core.transformation.collect.EqualsSplitCriteria;
import redberry.core.transformation.collect.ScalarsSplitCriteria;
import redberry.core.transformation.concurrent.ExpandAndCollectTransformation;
import redberry.core.transformation.contractions.IndexesContractionsTransformation;
import redberry.core.utils.Indicator;
import redberryphysics.core.util.IndexesFactoryUtil;
import static core.TAssert.*;
import static redberryphysics.core.util.IndexesFactoryUtil.*;

/**
 *
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 * @author Konstantin Kiselev
 */
public class OneLoop1Test {
    public OneLoop1Test() {
    }

    @Test
    public void evalRR() {
        OneLoop loop1 = new OneLoop(OneLoop.EVAL.INITIALIZE);
        loop1.evalRR();
        System.out.println(loop1.RR.toString(ToStringMode.UTF8));
    }

    @Test
    public void evalHATK() {
        OneLoop loop1 = new OneLoop(OneLoop.EVAL.INITIALIZE);
        loop1.evalHatK();
        System.out.println(loop1.HATK_1.toString(ToStringMode.UTF8));
    }

    @Test
    public void evalDeltas() {
        OneLoop loop1 = new OneLoop(OneLoop.EVAL.INITIALIZE);
        loop1.evalDeltas();
        Sum s = (Sum) loop1.DELTA_3.right();
        loop1.evalHatK();
        CollectInputPort cip = new CollectInputPort(EqualsSplitCriteria.INSTANCE1);
        System.out.println(s.size());
        Product p = (Product) s.getElements().get(6);
        for (Expression hatk : loop1.HATKs)
            hatk.asSubstitution().transform(p);

        Tensor p1 = new Product(p.getElements().get(1).clone(), p.getElements().get(2).clone());
        p1 = ExpandAndCollectTransformation.EXPAND_AND_COLLECT1.transform(p1);
        Tensor p2 = new Product(p1, p.getElements().get(3).clone());
        p2 = ExpandAndCollectTransformation.EXPAND_AND_COLLECT1.transform(p2);
        System.out.println(p2.toString(ToStringMode.UTF8));
//        for (int i = 0; i < s.size(); ++i) {
//            Tensor t = s.getElements().get(i);
//
//            for (Expression hatk : loop1.HATKs)
//                t = hatk.asSubstitution().transform(t);
//            OutputPortUnsafe<Tensor> bracketsOutput = ExpandBracketsOutput.create(t, Indicator.SYMBOL_INDICATOR);
//            Tensor cur;
//            while ((cur = bracketsOutput.take()) != null) {
//                cur = Transformations.contractMetrics(cur);
//                cur = loop1.KRONECKER_DIMENSION.asSubstitution().transform(cur);
//                cur = Transformations.calculateNumbers(cur);
//                cip.put(cur);
//            }
//            System.out.println(i);
//        }
    }

    @Ignore
    @Test
    public void test() {
        OneLoop loop1 = new OneLoop(OneLoop.EVAL.INITIALIZE);
        loop1.evalRR();
        assertIndexes(loop1.RR);
        loop1.evalDeltas();
        System.out.println(loop1.DELTA_3.toString(ToStringMode.UTF8));
        System.out.println(loop1.DELTA_2.toString(ToStringMode.UTF8));
        System.out.println(loop1.DELTA_1.toString(ToStringMode.UTF8));
        System.out.println(loop1.DELTA_4.toString(ToStringMode.UTF8));
        assertIndexes(loop1.DELTAs);
//        loop1.DELTA_3.asSubstitution();
        loop1.DELTA_4.asSubstitution();
    }

    @Ignore
    @Test
    public void testHATKs() {
        OneLoop loop1 = new OneLoop(OneLoop.EVAL.EVAL_HATK);
        System.out.println(loop1.HATK_1.toString(ToStringMode.UTF8));
        System.out.println(loop1.HATK_2.toString(ToStringMode.UTF8));
        System.out.println(loop1.DELTA_1.toString(ToStringMode.UTF8));
        System.out.println(loop1.DELTA_2.toString(ToStringMode.UTF8));
    }

    @Ignore
    @Test
    public void testAll() {
        OneLoop loop1 = new OneLoop(OneLoop.EVAL.EVAL_ALL);
        System.out.println(loop1.RR.toString(ToStringMode.UTF8));
        System.out.println(((MultiTensor) loop1.RR.right()).size());
        assertIndexes(loop1.RR);
    }

    @Ignore
    @Test
    public void testDelta2() {
        OneLoop loop1 = new OneLoop(OneLoop.EVAL.EVAL_HATK);
        IndexesInsertion indexesInsertion = new IndexesInsertion(loop1.matricesIndicator, IndexesFactoryUtil.createIndexes(loop1.DELTAs, "^{\\mu\\nu}_{\\alpha\\beta}"));

        loop1.DELTA_2.eval(
                indexesInsertion,
                loop1.L.asSubstitution(),
                CalculateNumbers.INSTANCE);

        System.out.println("HATK subs ... ");
        loop1.DELTA_2.eval(
                loop1.HATK_1.asSubstitution(),
                loop1.HATK_2.asSubstitution(),
                loop1.HATK_3.asSubstitution(),
                loop1.HATK_4.asSubstitution());
        System.out.println("HATK subs ... done");

        System.out.println("Expand brackets ...");
        loop1.DELTA_2.eval(
                new Transformer(ExpandBrackets.EXPAND_EXCEPT_SYMBOLS));
        System.out.println("Expand brackets ... done");

        System.out.println("Indexes contractions ...");
        loop1.DELTA_2.eval(
                IndexesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                loop1.KRONECKER_DIMENSION.asSubstitution(),
                CalculateNumbers.INSTANCE);
        System.out.println("Indexes contractions ... done");
        assertIndexes(loop1.DELTA_2.right());
        System.out.println("Collecting " + ((MultiTensor) loop1.DELTA_2.right()).size() + " elements ...");

        loop1.DELTA_2.eval(
                CollectFactory.createCollectEqualTerms());
        System.out.println("Collecting ... done. The new size: " + ((MultiTensor) loop1.DELTA_2.right()).size());

        System.out.println("Scalars ...");
        loop1.DELTA_2.eval(
                CalculateNumbers.INSTANCE,
                CollectFactory.createCollectAllScalars(),
                CalculateNumbers.INSTANCE);
        System.out.println("Scalars ... done");

        System.out.println(loop1.DELTA_2.toString(ToStringMode.UTF8));
        Tensor rhs = loop1.DELTA_2.right();
        TensorTreeIterator iterator = new TensorFirstTreeIterator(rhs);
        while (iterator.hasNext())
            if (TTest.testIsSymbol(iterator.next()))
                iterator.remove();
        System.out.println(rhs.toString(ToStringMode.UTF8));
    }

    @Ignore
    @Test
    public void test1() {
        OneLoop loop1 = new OneLoop(OneLoop.EVAL.INITIALIZE);
        loop1.evalRR();
        loop1.evalDeltas();
        for (Expression delta : loop1.DELTAs)
            loop1.RR.eval(delta.asSubstitution());
        loop1.RR.eval(
                IndexesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                loop1.KRONECKER_DIMENSION.asSubstitution(),
                CalculateNumbers.INSTANCE,
                new Transformer(ExpandBrackets.EXPAND_EXCEPT_SYMBOLS),
                IndexesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                CollectFactory.createCollectEqualTerms(),
                CollectFactory.createCollectScalars(),
                CalculateNumbers.INSTANCE);

        Transformation indexesInsertion;
        indexesInsertion = new IndexesInsertion(loop1.matricesIndicator, createIndexes(loop1.HATKs, "^{\\mu\\nu}_{\\alpha\\beta}"));
        for (Expression hatK : loop1.HATKs) {
            hatK.eval(indexesInsertion);
            loop1.RR.eval(hatK.asSubstitution());
        }


//        loop1.MATRIX_K.eval(
//                loop1.P.asSubstitution(),
//                new Transformer(ExpandBrackets.EXPAND_EXCEPT_SYMBOLS),
//                CollectFactory.createCollectEqualTerms(),
//                CollectFactory.createCollectEqualTerms());
        System.out.println("matrix K evaluating");

        loop1.MATRIX_K.eval(loop1.P.asSubstitution(),
                new Transformer(ExpandBrackets.EXPAND_EXCEPT_SYMBOLS),
                CollectFactory.createCollectEqualTerms(),
                CalculateNumbers.INSTANCE);
        loop1.MATRIX_K_INV.eval(loop1.P.asSubstitution(),
                new Transformer(ExpandBrackets.EXPAND_EXCEPT_SYMBOLS),
                CollectFactory.createCollectEqualTerms(),
                CalculateNumbers.INSTANCE);
        System.out.println("done");
        loop1.RR.eval(IndexesContractionsTransformation.CONTRACTIONS_WITH_METRIC);
        System.out.println(((MultiTensor) loop1.RR.right()).getElements().get(0).toString(ToStringMode.UTF8));
        loop1.RR.eval(loop1.MATRIX_K_INV.asSubstitution());//, CollectFactory.createCollectEqualTerms(), CalculateNumbers.INSTANCE, CollectFactory.createCollectAllScalars());
        loop1.RR.eval(loop1.MATRIX_K.asSubstitution());


        Tensor rhs = ((MultiTensor) loop1.RR.right()).getElements().get(0).clone();
        System.out.println(" Evaluating RR ");
        loop1.RR.eval(
                IndexesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                loop1.KRONECKER_DIMENSION.asSubstitution());

        TensorTreeIterator treeIterator = new TensorFirstTreeIterator(rhs);
        Tensor current;
        while (treeIterator.hasNext()) {
            current = treeIterator.next();
            if (current instanceof Sum && !TTest.testIsSymbol(current)) {
                current = CollectFactory.createCollectEqualTerms().transform(current.clone());
                current = CalculateNumbers.INSTANCE.transform(current);
                treeIterator.set(current);
            }
        }
        List<Sum> sums = new ArrayList<>();
        TensorIterator iterator = rhs.iterator();
        while (iterator.hasNext()) {
            current = iterator.next();
            if (current instanceof Sum && !TTest.testIsSymbol(current)) {
                sums.add((Sum) current.clone());
                iterator.remove();
            }
        }


        System.out.println(" counting terms: ... ");
        long start = System.currentTimeMillis();
        Transformation ec = new ExpandAndCollectTransformation(
                Indicator.SYMBOL_INDICATOR,
                new Transformation[]{IndexesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                    loop1.KRONECKER_DIMENSION.asSubstitution(),
                    CalculateNumbers.INSTANCE});
        for (Sum sum : sums) {
            rhs = new Product(rhs.equivalent(), sum);
            rhs = ec.transform(rhs.clone());
            System.out.println("+");
        }
        rhs = loop1.KRONECKER_DIMENSION.asSubstitution().transform(rhs);
        rhs = Transformations.calculateNumbers(rhs);
        long stop = System.currentTimeMillis();
        System.out.println("... . Evalution time = " + (stop - start) + "ms");
        System.out.println(rhs.toString(ToStringMode.UTF8));
//        System.out.println(fS.transform(rhs).toString(ToStringMode.UTF8));

//        Tensor rhs_1 = ((MultiTensor) rhs).getElements().get(0);
//
//        System.out.println("contracting  ");
//        rhs_1 = IndexesContractionsTransformation.CONTRACTIONS_WITH_METRIC.transform(rhs_1);
//        rhs_1 = loop1.KRONECKER_DIMENSION.asSubstitution().transform(rhs_1);
//        rhs_1 = CalculateNumbers.INSTANCE.transform(rhs_1);
//        System.out.println("done");
//        System.out.println("expanding  ");
//        rhs_1 = new Transformer(ExpandBrackets.EXPAND_EXCEPT_SYMBOLS).transform(rhs_1);
//        System.out.println("done");
//        System.out.println("contracting  ");
//        rhs_1 = IndexesContractionsTransformation.CONTRACTIONS_WITH_METRIC.transform(rhs_1);
//        rhs_1 = loop1.KRONECKER_DIMENSION.asSubstitution().transform(rhs_1);
//        rhs_1 = CalculateNumbers.INSTANCE.transform(rhs_1);
//        System.out.println("done");
//        System.out.println(rhs_1.toString(ToStringMode.UTF8));
//
//        int oldSize = ((MultiTensor) rhs_1).size();
//        TensorIterator iterator = rhs_1.iterator();
//        Tensor c;
//        int count;
//        while (iterator.hasNext()) {
//            c = iterator.next();
//            count = 0;
//            for (Tensor n : c)
//                if (TTest.testEqualstensorStructure(n, CC.parse("n_{\\mu}")))
//                    count++;
//            if (count % 2 != 0)
//                iterator.remove();
//        }
//        int newSize = ((MultiTensor) rhs_1).size();
//        System.out.println("Old size: " + oldSize + ", new sizeL " + newSize + ".  Collecting");
//        Tensor nT = ((MultiTensor) rhs_1).getElements().get(0);
////        count = 0;
////        for (Tensor t : rhs_1) {
////            if (count > 10)
////                break;
////            ((Sum) nT).add(t);
////        }
//        System.out.println(nT.toString(ToStringMode.UTF8));
//        nT = loop1.MATRIX_K.asSubstitution().transform(nT);
//        nT = Transformations.expandBrackets(nT);
//        nT = Transformations.contractMetrics(nT);
//        nT = loop1.KRONECKER_DIMENSION.asSubstitution().transform(nT);
//        nT = Transformations.calculateNumbers(nT);
//        System.out.println("Collecting");
//
//        nT = CollectFactory.createCollectEqualTerms().transform(nT);
//        System.out.println(((Sum) nT).size());
//        System.out.println(nT.toString(ToStringMode.UTF8));
    }

//    @Ignore
    @Test
    public void test2() {
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
//        System.out.println(rhs.toString(ToStringMode.UTF8));
        System.out.println(" Evaluating RR ");

        TensorTreeIterator treeIterator = new TensorFirstTreeIterator(rhs);
        Tensor current;
//        while (treeIterator.hasNext()) {
//            current = treeIterator.next();
//            if (current instanceof Sum && !TTest.testIsSymbol(current)) {
//                current = CollectFactory.createCollectEqualTerms().transform(current.clone());
//                current = CalculateNumbers.INSTANCE.transform(current);
//                treeIterator.set(current);
//            }
//        }
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
        Product firstPart = new Product(sums.get(1), sums.get(2));
        System.out.println(sums.get(1).size());
        System.out.println(sums.get(2).size());
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

    @Ignore
    @Test
    public void test2_2() {
        OneLoop loop1 = new OneLoop(OneLoop.EVAL.INITIALIZE);
        loop1.evalRR();
        loop1.evalHatK();
        loop1.evalDeltas();
        for (Expression delta : loop1.DELTAs)
            loop1.RR.eval(delta.asSubstitution());
        System.out.println(loop1.RR.toString(ToStringMode.UTF8));
        for (Expression hatK : loop1.HATKs)
            loop1.RR.eval(hatK.asSubstitution());
        loop1.RR.eval(
                IndexesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                loop1.KRONECKER_DIMENSION.asSubstitution());

        int sC = 0, tC = 0;;
        for (Tensor t : loop1.RR.right()) {
            sC = 0;
            for (Tensor m : t)
                if (m instanceof Sum && !TTest.testIsSymbol(m))
                    sC++;
//            if (sC != 4)
            System.out.println(sC);
            tC++;
        }
        System.out.println(tC);
    }
}
