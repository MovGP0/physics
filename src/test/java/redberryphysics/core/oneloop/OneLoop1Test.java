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

import com.sun.org.apache.bcel.internal.generic.FASTORE;
import java.util.ArrayList;
import java.util.List;
import org.junit.Ignore;
import redberry.core.tensor.Sum;
import redberry.core.tensor.Expression;
import redberry.core.transformation.Transformation;
import redberry.core.tensor.Tensor;
import org.junit.Test;
import redberry.concurrent.OutputPortUnsafe;
import redberry.core.context.ToStringMode;
import redberry.core.tensor.MultiTensor;
import redberry.core.tensor.Product;
import redberry.core.tensor.TensorIterator;
import redberry.core.tensor.iterators.TensorFirstTreeIterator;
import redberry.core.tensor.iterators.TensorTreeIterator;
import redberry.core.tensor.test.TTest;
import redberry.core.transformation.CalculateNumbers;
import redberry.core.transformation.ExpandBrackets;
import redberry.core.transformation.IndexesInsertion;
import redberry.core.transformation.Transformations;
import redberry.core.transformation.Transformer;
import redberry.core.transformation.collect.AbstractCollectTerms;
import redberry.core.transformation.collect.CollectFactory;
import redberry.core.transformation.collect.EqualsSplitCriteria;
import redberry.core.transformation.concurrent.ExpandBracketsOutput;
import redberry.core.transformation.concurrent.ExpandBracketsOutputTransformation;
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

    @Ignore
    @Test
    public void test() {
        OneLoop1 loop1 = new OneLoop1(OneLoop1.EVAL.INITIALIZE);
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
        OneLoop1 loop1 = new OneLoop1(OneLoop1.EVAL.EVAL_HATK);
        System.out.println(loop1.HATK_1.toString(ToStringMode.UTF8));
        System.out.println(loop1.HATK_2.toString(ToStringMode.UTF8));
        System.out.println(loop1.DELTA_1.toString(ToStringMode.UTF8));
        System.out.println(loop1.DELTA_2.toString(ToStringMode.UTF8));
    }

    @Ignore
    @Test
    public void testAll() {
        OneLoop1 loop1 = new OneLoop1(OneLoop1.EVAL.EVAL_ALL);
        System.out.println(loop1.RR.toString(ToStringMode.UTF8));
        System.out.println(((MultiTensor) loop1.RR.right()).size());
        assertIndexes(loop1.RR);
    }

    @Ignore
    @Test
    public void testDelta2() {
        OneLoop1 loop1 = new OneLoop1(OneLoop1.EVAL.EVAL_HATK);
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
                CollectFactory.ccreateCollectAllScalars(),
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

    @Test
    public void test1() {
        OneLoop1 loop1 = new OneLoop1(OneLoop1.EVAL.INITIALIZE);
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
//        loop1.MATRIX_K.eval(loop1.P.asSubstitution(), CollectFactory.createCollectEqualTerms(), CalculateNumbers.INSTANCE);
//        loop1.MATRIX_K_INV.eval(loop1.P.asSubstitution(), CollectFactory.createCollectEqualTerms(), CalculateNumbers.INSTANCE);
        System.out.println("done");
        loop1.RR.eval(loop1.MATRIX_K_INV.asSubstitution());//, CollectFactory.createCollectEqualTerms(), CalculateNumbers.INSTANCE, CollectFactory.ccreateCollectAllScalars());
        loop1.RR.eval(loop1.MATRIX_K.asSubstitution());


        Tensor rhs = ((MultiTensor) loop1.RR.right()).getElements().get(0).clone();
        System.out.println(" Evaluating RR ");

        List<Sum> sums = new ArrayList<>();
        TensorIterator iterator = rhs.iterator();
        Tensor current;
        while (iterator.hasNext()) {
            current = iterator.next();
            if (current instanceof Sum && !TTest.testIsSymbol(current)) {
                sums.add((Sum) current.clone());
                iterator.remove();
            }
        }


        System.out.println(" counting terms: ... ");
        long start = System.currentTimeMillis();
        ExpandBracketsOutputTransformation fS = ExpandBracketsOutputTransformation.EXPAND_COMMON;
        for (Sum sum : sums) {
            rhs = new Product(rhs.equivalent(), sum);
            rhs = fS.transform(rhs.clone());
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

    @Ignore
    @Test
    public void test2() {
        OneLoop1 loop1 = new OneLoop1(OneLoop1.EVAL.INITIALIZE);
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
        System.out.println("Evaluating \\hat K ...");
        loop1.evalHatK();
        System.out.println("Done: \\hat K_1 size =" + ((MultiTensor) loop1.HATK_1.right()).size() + "  \\hat K_2 size = " + ((MultiTensor) loop1.HATK_2.right()).size());
        Tensor term = ((MultiTensor) loop1.RR.right()).getElements().get(0);
        System.out.println(term.toString(ToStringMode.UTF8));
        term = loop1.HATK_1.asSubstitution().transform(term);
//        term = loop1.HATK_2.asSubstitution().transform(term);
        term = loop1.HATK_3.asSubstitution().transform(term);
        System.out.println("Expanding...");
        term = new Transformer(ExpandBrackets.EXPAND_EXCEPT_SYMBOLS).transform(term);
        System.out.println("Contracting...");
        term = Transformations.contractMetrics(term);
        System.out.println("Done. Result sum size = " + ((MultiTensor) term).size());

        term = ((Sum) term).getElements().get(0);
//        term = loop1.HATK_2.asSubstitution().transform(term);
        term = loop1.HATK_2.asSubstitution().transform(term);
        System.out.println("Expanding...");
        term = new Transformer(ExpandBrackets.EXPAND_EXCEPT_SYMBOLS).transform(term);
        System.out.println("Contracting...");
        term = Transformations.contractMetrics(term);
        System.out.println("Done. Result sum size = " + ((MultiTensor) term).size());
        System.out.println("Collecting...");
        term = CollectFactory.createCollectEqualTerms().transform(term);
        System.out.println("Collecting scalars...");
        term = CollectFactory.createCollectScalars().transform(term);

    }
}
