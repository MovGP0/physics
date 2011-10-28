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

import redberry.core.tensor.Sum;
import redberry.core.tensor.Tensor;
import org.junit.BeforeClass;
import org.junit.Test;
import redberry.core.context.CC;
import redberry.core.context.ToStringMode;
import redberry.core.tensor.MultiTensor;
import redberry.core.tensor.iterators.TensorFirstTreeIterator;
import redberry.core.tensor.iterators.TensorTreeIterator;
import redberry.core.tensor.test.TTest;
import redberry.core.transformation.CalculateNumbers;
import redberry.core.transformation.ExpandBrackets;
import redberry.core.transformation.IndexesInsertion;
import redberry.core.transformation.Transformations;
import redberry.core.transformation.Transformer;
import redberry.core.transformation.collect.CollectFactory;
import redberry.core.transformation.contractions.IndexesContractionsTransformation;
import redberry.core.utils.TensorUtils;
import redberryphysics.core.util.IndexesFactoryUtil;
import static core.TAssert.*;
import static org.junit.Assert.*;
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

    @Test
    public void testHATKs() {
        OneLoop1 loop1 = new OneLoop1(OneLoop1.EVAL.EVAL_HATK);
        System.out.println(loop1.HATK_1.toString(ToStringMode.UTF8));
        System.out.println(loop1.HATK_2.toString(ToStringMode.UTF8));
        System.out.println(loop1.DELTA_1.toString(ToStringMode.UTF8));
        System.out.println(loop1.DELTA_2.toString(ToStringMode.UTF8));
    }

    @Test
    public void testAll() {
        OneLoop1 loop1 = new OneLoop1(OneLoop1.EVAL.EVAL_ALL);
        System.out.println(loop1.RR.toString(ToStringMode.UTF8));
        System.out.println(((MultiTensor) loop1.RR.right()).size());
//        Tensor t = loop1.RR.right();
//        for (Tensor f : t)
//            if (!TensorUtils.testIndexesConsistent(f))
//                System.out.println(f);
//        TensorTreeIterator iterator  = new TensorFirstTreeIterator(t);
//        Tensor next;
//        while(iterator.hasNext())
//        {
//            next = iterator.next();
//            iterator.set(Transformations.renameConflictingIndexes(next));
//        }
        assertIndexes(loop1.RR);
    }

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
        Sum t = (Sum) CC.parse("g^{μν}*d_{ε}^{γ}*d_{ζ}^{δ}+g^{μν}*d_{ε}^{δ}*d_{ζ}^{γ}+g^{γδ}*g^{μν}*g_{ζε}+n_{ε}*n^{γ}*g^{μν}*d_{ζ}^{δ}+n_{ε}*n^{δ}*d_{ζ}^{γ}*g^{μν}+n_{ζ}*n^{γ}*g^{μν}*d_{ε}^{δ}+n_{ζ}*n^{δ}*d_{ε}^{γ}*g^{μν}+n_{ε}*n_{ζ}*g^{γδ}*g^{μν}+n^{γ}*n^{δ}*g^{μν}*g_{ζε}+n_{ε}*n_{ζ}*n^{γ}*n^{δ}*g^{μν}+n^{β}*n_{β}*g^{γδ}*g^{μν}*g_{εζ}+n^{β}*n_{β}*n^{γ}*n^{δ}*g^{μν}*g_{εζ}+g^{νγ}*d_{ζ}^{δ}*d^{μ}_{ε}+d_{ζ}^{γ}*g^{νδ}*d^{μ}_{ε}+g^{γδ}*d^{μ}_{ε}*d_{ζ}^{ν}+n^{ν}*n^{γ}*d_{ζ}^{δ}*d^{μ}_{ε}+n^{ν}*n^{δ}*d_{ζ}^{γ}*d^{μ}_{ε}+n_{ζ}*n^{γ}*g^{νδ}*d^{μ}_{ε}+n_{ζ}*n^{δ}*g^{νγ}*d^{μ}_{ε}+n^{ν}*n_{ζ}*g^{γδ}*d^{μ}_{ε}+n^{γ}*n^{δ}*d^{μ}_{ε}*d_{ζ}^{ν}+n^{ν}*n_{ζ}*n^{γ}*n^{δ}*d^{μ}_{ε}+g^{νγ}*d_{ε}^{δ}*d^{μ}_{ζ}+d_{ε}^{γ}*g^{νδ}*d^{μ}_{ζ}+g^{γδ}*d^{μ}_{ζ}*d_{ε}^{ν}+n^{ν}*n^{γ}*d_{ε}^{δ}*d^{μ}_{ζ}+n^{ν}*n^{δ}*d_{ε}^{γ}*d^{μ}_{ζ}+n_{ε}*n^{γ}*g^{νδ}*d^{μ}_{ζ}+n_{ε}*n^{δ}*g^{νγ}*d^{μ}_{ζ}+n^{ν}*n_{ε}*g^{γδ}*d^{μ}_{ζ}+n^{γ}*n^{δ}*d^{μ}_{ζ}*d_{ε}^{ν}+n^{ν}*n_{ε}*n^{γ}*n^{δ}*d^{μ}_{ζ}+g^{μγ}*d_{ζ}^{δ}*d^{ν}_{ε}+d_{ζ}^{γ}*g^{μδ}*d^{ν}_{ε}+n^{μ}*n^{γ}*d_{ζ}^{δ}*d^{ν}_{ε}+n^{μ}*n^{δ}*d_{ζ}^{γ}*d^{ν}_{ε}+n_{ζ}*n^{γ}*g^{μδ}*d^{ν}_{ε}+n_{ζ}*n^{δ}*g^{μγ}*d^{ν}_{ε}+n^{μ}*n_{ζ}*g^{γδ}*d^{ν}_{ε}+n^{μ}*n_{ζ}*n^{γ}*n^{δ}*d^{ν}_{ε}+g^{μγ}*d_{ε}^{δ}*d^{ν}_{ζ}+d_{ε}^{γ}*g^{μδ}*d^{ν}_{ζ}+n^{μ}*n^{γ}*d_{ε}^{δ}*d^{ν}_{ζ}+n^{μ}*n^{δ}*d_{ε}^{γ}*d^{ν}_{ζ}+n_{ε}*n^{γ}*g^{μδ}*d^{ν}_{ζ}+n_{ε}*n^{δ}*g^{μγ}*d^{ν}_{ζ}+n^{μ}*n_{ε}*g^{γδ}*d^{ν}_{ζ}+n^{μ}*n_{ε}*n^{γ}*n^{δ}*d^{ν}_{ζ}+g^{μγ}*g^{νδ}*g_{εζ}+g^{νγ}*g^{μδ}*g_{εζ}+n^{μ}*n^{γ}*g^{νδ}*g_{εζ}+n^{μ}*n^{δ}*g^{νγ}*g_{εζ}+n^{ν}*n^{γ}*g^{μδ}*g_{εζ}+n^{ν}*n^{δ}*g^{μγ}*g_{εζ}+n^{μ}*n^{ν}*g^{γδ}*g_{εζ}+n^{μ}*n^{ν}*n^{γ}*n^{δ}*g_{εζ}+n^{β}*n_{β}*g^{γδ}*d^{μ}_{ε}*d^{ν}_{ζ}+n^{β}*n_{β}*n^{γ}*n^{δ}*d^{μ}_{ε}*d^{ν}_{ζ}+n^{β}*n_{β}*g^{γδ}*d^{μ}_{ζ}*d^{ν}_{ε}+n^{β}*n_{β}*n^{γ}*n^{δ}*d^{μ}_{ζ}*d^{ν}_{ε}+n^{μ}*n^{ν}*d_{ε}^{γ}*d_{ζ}^{δ}+n^{μ}*n^{ν}*d_{ε}^{δ}*d_{ζ}^{γ}+n_{ε}*n^{γ}*n^{μ}*n^{ν}*d_{ζ}^{δ}+n_{ε}*n^{δ}*n^{μ}*n^{ν}*d_{ζ}^{γ}+n_{ζ}*n^{γ}*n^{μ}*n^{ν}*d_{ε}^{δ}+n_{ζ}*n^{δ}*n^{μ}*n^{ν}*d_{ε}^{γ}+n_{ε}*n_{ζ}*n^{μ}*n^{ν}*g^{γδ}+n_{ε}*n_{ζ}*n^{γ}*n^{δ}*n^{μ}*n^{ν}+n^{ι}*n_{ι}*n^{μ}*n^{ν}*g^{γδ}*g_{ζε}+n^{ι}*n_{ι}*n^{γ}*n^{δ}*n^{μ}*n^{ν}*g_{ζε}+n^{λ}*n^{γ}*n_{λ}*n^{ν}*d_{ζ}^{δ}*d_{ε}^{μ}+n^{λ}*n^{δ}*n_{λ}*n^{ν}*d_{ζ}^{γ}*d_{ε}^{μ}+n^{λ}*n_{ζ}*n_{λ}*n^{ν}*g^{γδ}*d_{ε}^{μ}+n^{λ}*n_{ζ}*n^{γ}*n^{δ}*n_{λ}*n^{ν}*d_{ε}^{μ}+n^{λ}*n^{γ}*n_{λ}*n^{ν}*d_{ε}^{δ}*d_{ζ}^{μ}+n^{λ}*n^{δ}*n_{λ}*n^{ν}*d_{ε}^{γ}*d_{ζ}^{μ}+n^{λ}*n_{ε}*n_{λ}*n^{ν}*g^{γδ}*d_{ζ}^{μ}+n^{λ}*n_{ε}*n^{γ}*n^{δ}*n_{λ}*n^{ν}*d_{ζ}^{μ}+n_{ε}*n^{ν}*g^{μγ}*d_{ζ}^{δ}+n_{ε}*n^{ν}*d_{ζ}^{γ}*g^{μδ}+n_{ζ}*n^{γ}*n_{ε}*n^{ν}*g^{μδ}+n_{ζ}*n^{δ}*n_{ε}*n^{ν}*g^{μγ}+n_{ζ}*n^{ν}*g^{μγ}*d_{ε}^{δ}+n_{ζ}*n^{ν}*d_{ε}^{γ}*g^{μδ}+n^{λ}*n^{γ}*n_{λ}*n^{ν}*g^{μδ}*g_{ζε}+n^{λ}*n^{δ}*n_{λ}*n^{ν}*g^{μγ}*g_{ζε}+n_{η}*n^{γ}*n^{μ}*n_{ε}*n^{η}*n^{ν}*d_{ζ}^{δ}+n_{η}*n^{δ}*n^{μ}*n_{ε}*n^{η}*n^{ν}*d_{ζ}^{γ}+n_{η}*n_{ζ}*n^{μ}*n_{ε}*n^{η}*n^{ν}*g^{γδ}+n_{η}*n_{ζ}*n^{γ}*n^{δ}*n^{μ}*n_{ε}*n^{η}*n^{ν}+n^{λ}*n_{η}*n_{λ}*n_{ε}*n^{η}*n^{ν}*g^{γδ}*d_{ζ}^{μ}+n^{λ}*n_{η}*n^{γ}*n^{δ}*n_{λ}*n_{ε}*n^{η}*n^{ν}*d_{ζ}^{μ}+n_{η}*n_{ε}*n^{η}*n^{ν}*g^{μγ}*d_{ζ}^{δ}+n_{η}*n_{ε}*n^{η}*n^{ν}*d_{ζ}^{γ}*g^{μδ}+n_{ζ}*n^{γ}*n_{η}*n_{ε}*n^{η}*n^{ν}*g^{μδ}+n_{ζ}*n^{δ}*n_{η}*n_{ε}*n^{η}*n^{ν}*g^{μγ}+n_{η}*n^{γ}*n^{μ}*n_{ζ}*n^{η}*n^{ν}*d_{ε}^{δ}+n_{η}*n^{δ}*n^{μ}*n_{ζ}*n^{η}*n^{ν}*d_{ε}^{γ}+n^{λ}*n_{η}*n_{λ}*n_{ζ}*n^{η}*n^{ν}*g^{γδ}*d_{ε}^{μ}+n^{λ}*n_{η}*n^{γ}*n^{δ}*n_{λ}*n_{ζ}*n^{η}*n^{ν}*d_{ε}^{μ}+n_{η}*n_{ζ}*n^{η}*n^{ν}*g^{μγ}*d_{ε}^{δ}+n_{η}*n_{ζ}*n^{η}*n^{ν}*d_{ε}^{γ}*g^{μδ}+n_{η}*n_{θ}*n^{μ}*n^{η}*n^{θ}*n^{ν}*g^{γδ}*g_{ζε}+n_{η}*n_{θ}*n^{γ}*n^{δ}*n^{μ}*n^{η}*n^{θ}*n^{ν}*g_{ζε}+n_{θ}*n^{γ}*n_{η}*n^{η}*n^{θ}*n^{ν}*g^{μδ}*g_{ζε}+n_{θ}*n^{δ}*n_{η}*n^{η}*n^{θ}*n^{ν}*g^{μγ}*g_{ζε}+n_{η}*n_{θ}*n^{μ}*n_{ε}*n_{ζ}*n^{η}*n^{θ}*n^{ν}*g^{γδ}+n_{η}*n_{θ}*n^{γ}*n^{δ}*n^{μ}*n_{ε}*n_{ζ}*n^{η}*n^{θ}*n^{ν}+n_{θ}*n^{γ}*n_{η}*n_{ε}*n_{ζ}*n^{η}*n^{θ}*n^{ν}*g^{μδ}+n_{θ}*n^{δ}*n_{η}*n_{ε}*n_{ζ}*n^{η}*n^{θ}*n^{ν}*g^{μγ}+n_{η}*n_{θ}*n^{μ}*n^{ξ}*n_{ξ}*n^{η}*n^{θ}*n^{ν}*g^{γδ}*g_{εζ}+n_{η}*n_{θ}*n^{γ}*n^{δ}*n^{μ}*n^{ξ}*n_{ξ}*n^{η}*n^{θ}*n^{ν}*g_{εζ}+n_{θ}*n^{γ}*n_{η}*n^{ξ}*n_{ξ}*n^{η}*n^{θ}*n^{ν}*g^{μδ}*g_{εζ}+n_{θ}*n^{δ}*n_{η}*n^{ξ}*n_{ξ}*n^{η}*n^{θ}*n^{ν}*g^{μγ}*g_{εζ}+n_{η}*n^{γ}*n^{μ}*n^{η}*d_{ζ}^{δ}*d^{ν}_{ε}+n_{η}*n^{δ}*n^{μ}*n^{η}*d_{ζ}^{γ}*d^{ν}_{ε}+n_{η}*n_{ζ}*n^{μ}*n^{η}*g^{γδ}*d^{ν}_{ε}+n_{η}*n_{ζ}*n^{γ}*n^{δ}*n^{μ}*n^{η}*d^{ν}_{ε}+n^{λ}*n_{η}*n_{λ}*n^{η}*g^{γδ}*d_{ζ}^{μ}*d^{ν}_{ε}+n^{λ}*n_{η}*n^{γ}*n^{δ}*n_{λ}*n^{η}*d_{ζ}^{μ}*d^{ν}_{ε}+n_{η}*n^{η}*g^{μγ}*d_{ζ}^{δ}*d^{ν}_{ε}+n_{η}*n^{η}*d_{ζ}^{γ}*g^{μδ}*d^{ν}_{ε}+n_{ζ}*n^{γ}*n_{η}*n^{η}*g^{μδ}*d^{ν}_{ε}+n_{ζ}*n^{δ}*n_{η}*n^{η}*g^{μγ}*d^{ν}_{ε}+n_{η}*n^{γ}*n^{μ}*n^{α}*n^{η}*n_{α}*d_{ζ}^{δ}*d^{ν}_{ε}+n_{η}*n^{δ}*n^{μ}*n^{α}*n^{η}*n_{α}*d_{ζ}^{γ}*d^{ν}_{ε}+n_{η}*n_{ζ}*n^{μ}*n^{α}*n^{η}*n_{α}*g^{γδ}*d^{ν}_{ε}+n_{η}*n_{ζ}*n^{γ}*n^{δ}*n^{μ}*n^{α}*n^{η}*n_{α}*d^{ν}_{ε}+n^{λ}*n_{η}*n_{λ}*n^{α}*n^{η}*n_{α}*g^{γδ}*d_{ζ}^{μ}*d^{ν}_{ε}+n^{λ}*n_{η}*n^{γ}*n^{δ}*n_{λ}*n^{α}*n^{η}*n_{α}*d_{ζ}^{μ}*d^{ν}_{ε}+n_{η}*n^{α}*n^{η}*n_{α}*g^{μγ}*d_{ζ}^{δ}*d^{ν}_{ε}+n_{η}*n^{α}*n^{η}*n_{α}*d_{ζ}^{γ}*g^{μδ}*d^{ν}_{ε}+n_{ζ}*n^{γ}*n_{η}*n^{α}*n^{η}*n_{α}*g^{μδ}*d^{ν}_{ε}+n_{ζ}*n^{δ}*n_{η}*n^{α}*n^{η}*n_{α}*g^{μγ}*d^{ν}_{ε}+n_{η}*n_{θ}*n^{μ}*n^{α}*n_{ζ}*n^{η}*n^{θ}*n_{α}*g^{γδ}*d^{ν}_{ε}+n_{η}*n_{θ}*n^{γ}*n^{δ}*n^{μ}*n^{α}*n_{ζ}*n^{η}*n^{θ}*n_{α}*d^{ν}_{ε}+n_{θ}*n^{γ}*n_{η}*n^{α}*n_{ζ}*n^{η}*n^{θ}*n_{α}*g^{μδ}*d^{ν}_{ε}+n_{θ}*n^{δ}*n_{η}*n^{α}*n_{ζ}*n^{η}*n^{θ}*n_{α}*g^{μγ}*d^{ν}_{ε}+n_{η}*n^{γ}*n^{μ}*n^{η}*d_{ε}^{δ}*d^{ν}_{ζ}+n_{η}*n^{δ}*n^{μ}*n^{η}*d_{ε}^{γ}*d^{ν}_{ζ}+n_{η}*n_{ε}*n^{μ}*n^{η}*g^{γδ}*d^{ν}_{ζ}+n_{η}*n_{ε}*n^{γ}*n^{δ}*n^{μ}*n^{η}*d^{ν}_{ζ}+n^{λ}*n_{η}*n_{λ}*n^{η}*g^{γδ}*d_{ε}^{μ}*d^{ν}_{ζ}+n^{λ}*n_{η}*n^{γ}*n^{δ}*n_{λ}*n^{η}*d_{ε}^{μ}*d^{ν}_{ζ}+n_{η}*n^{η}*g^{μγ}*d_{ε}^{δ}*d^{ν}_{ζ}+n_{η}*n^{η}*d_{ε}^{γ}*g^{μδ}*d^{ν}_{ζ}+n_{ε}*n^{γ}*n_{η}*n^{η}*g^{μδ}*d^{ν}_{ζ}+n_{ε}*n^{δ}*n_{η}*n^{η}*g^{μγ}*d^{ν}_{ζ}+n_{η}*n^{γ}*n^{μ}*n^{α}*n^{η}*n_{α}*d_{ε}^{δ}*d^{ν}_{ζ}+n_{η}*n^{δ}*n^{μ}*n^{α}*n^{η}*n_{α}*d_{ε}^{γ}*d^{ν}_{ζ}+n_{η}*n_{ε}*n^{μ}*n^{α}*n^{η}*n_{α}*g^{γδ}*d^{ν}_{ζ}+n_{η}*n_{ε}*n^{γ}*n^{δ}*n^{μ}*n^{α}*n^{η}*n_{α}*d^{ν}_{ζ}+n^{λ}*n_{η}*n_{λ}*n^{α}*n^{η}*n_{α}*g^{γδ}*d_{ε}^{μ}*d^{ν}_{ζ}+n^{λ}*n_{η}*n^{γ}*n^{δ}*n_{λ}*n^{α}*n^{η}*n_{α}*d_{ε}^{μ}*d^{ν}_{ζ}+n_{η}*n^{α}*n^{η}*n_{α}*g^{μγ}*d_{ε}^{δ}*d^{ν}_{ζ}+n_{η}*n^{α}*n^{η}*n_{α}*d_{ε}^{γ}*g^{μδ}*d^{ν}_{ζ}+n_{ε}*n^{γ}*n_{η}*n^{α}*n^{η}*n_{α}*g^{μδ}*d^{ν}_{ζ}+n_{ε}*n^{δ}*n_{η}*n^{α}*n^{η}*n_{α}*g^{μγ}*d^{ν}_{ζ}+n_{η}*n_{θ}*n^{μ}*n^{α}*n_{ε}*n^{η}*n^{θ}*n_{α}*g^{γδ}*d^{ν}_{ζ}+n_{η}*n_{θ}*n^{γ}*n^{δ}*n^{μ}*n^{α}*n_{ε}*n^{η}*n^{θ}*n_{α}*d^{ν}_{ζ}+n_{θ}*n^{γ}*n_{η}*n^{α}*n_{ε}*n^{η}*n^{θ}*n_{α}*g^{μδ}*d^{ν}_{ζ}+n_{θ}*n^{δ}*n_{η}*n^{α}*n_{ε}*n^{η}*n^{θ}*n_{α}*g^{μγ}*d^{ν}_{ζ}+n^{μ}*n_{ε}*g^{νγ}*d_{ζ}^{δ}+n^{μ}*n_{ε}*d_{ζ}^{γ}*g^{νδ}+n_{ζ}*n^{γ}*n^{μ}*n_{ε}*g^{νδ}+n_{ζ}*n^{δ}*n^{μ}*n_{ε}*g^{νγ}+n^{λ}*n^{γ}*n_{λ}*n_{ε}*d_{ζ}^{δ}*g^{νμ}+n^{λ}*n^{δ}*n_{λ}*n_{ε}*d_{ζ}^{γ}*g^{νμ}+n^{λ}*n_{ζ}*n_{λ}*n_{ε}*g^{γδ}*g^{νμ}+n^{λ}*n_{ζ}*n^{γ}*n^{δ}*n_{λ}*n_{ε}*g^{νμ}+n^{λ}*n^{γ}*n_{λ}*n_{ε}*g^{νδ}*d_{ζ}^{μ}+n^{λ}*n^{δ}*n_{λ}*n_{ε}*g^{νγ}*d_{ζ}^{μ}+n_{ζ}*n_{ε}*g^{μγ}*g^{νδ}+n_{ζ}*n_{ε}*g^{νγ}*g^{μδ}+n_{η}*n^{γ}*n^{μ}*n_{ζ}*n^{η}*n_{ε}*g^{νδ}+n_{η}*n^{δ}*n^{μ}*n_{ζ}*n^{η}*n_{ε}*g^{νγ}+n^{λ}*n_{η}*n_{λ}*n_{ζ}*n^{η}*n_{ε}*g^{γδ}*g^{νμ}+n^{λ}*n_{η}*n^{γ}*n^{δ}*n_{λ}*n_{ζ}*n^{η}*n_{ε}*g^{νμ}+n_{η}*n_{ζ}*n^{η}*n_{ε}*g^{μγ}*g^{νδ}+n_{η}*n_{ζ}*n^{η}*n_{ε}*g^{νγ}*g^{μδ}+n^{μ}*n_{ζ}*g^{νγ}*d_{ε}^{δ}+n^{μ}*n_{ζ}*d_{ε}^{γ}*g^{νδ}+n^{λ}*n^{γ}*n_{λ}*n_{ζ}*d_{ε}^{δ}*g^{νμ}+n^{λ}*n^{δ}*n_{λ}*n_{ζ}*d_{ε}^{γ}*g^{νμ}+n^{λ}*n^{γ}*n_{λ}*n_{ζ}*g^{νδ}*d_{ε}^{μ}+n^{λ}*n^{δ}*n_{λ}*n_{ζ}*g^{νγ}*d_{ε}^{μ}+n_{θ}*n^{γ}*n^{μ}*n^{θ}*g^{νδ}*g_{εζ}+n_{θ}*n^{δ}*n^{μ}*n^{θ}*g^{νγ}*g_{εζ}+n^{λ}*n_{θ}*n_{λ}*n^{θ}*g^{γδ}*g^{νμ}*g_{εζ}+n^{λ}*n_{θ}*n^{γ}*n^{δ}*n_{λ}*n^{θ}*g^{νμ}*g_{εζ}+n_{θ}*n^{θ}*g^{μγ}*g^{νδ}*g_{εζ}+n_{θ}*n^{θ}*g^{νγ}*g^{μδ}*g_{εζ}+n_{η}*n^{γ}*n^{μ}*n^{α}*n^{η}*n_{α}*g^{νδ}*g_{εζ}+n_{η}*n^{δ}*n^{μ}*n^{α}*n^{η}*n_{α}*g^{νγ}*g_{εζ}+n^{λ}*n_{η}*n_{λ}*n^{α}*n^{η}*n_{α}*g^{γδ}*g^{νμ}*g_{εζ}+n^{λ}*n_{η}*n^{γ}*n^{δ}*n_{λ}*n^{α}*n^{η}*n_{α}*g^{νμ}*g_{εζ}+n_{η}*n^{α}*n^{η}*n_{α}*g^{μγ}*g^{νδ}*g_{εζ}+n_{η}*n^{α}*n^{η}*n_{α}*g^{νγ}*g^{μδ}*g_{εζ}+n_{λ}*n_{ε}*n^{λ}*n^{μ}*g^{νγ}*d_{ζ}^{δ}+n_{λ}*n_{ε}*n^{λ}*n^{μ}*d_{ζ}^{γ}*g^{νδ}+n_{λ}*n_{ζ}*n^{λ}*n^{μ}*g^{νγ}*d_{ε}^{δ}+n_{λ}*n_{ζ}*n^{λ}*n^{μ}*d_{ε}^{γ}*g^{νδ}+n_{ξ}*n^{γ}*n_{λ}*n_{ε}*n_{ζ}*n^{λ}*n^{ξ}*n^{μ}*g^{νδ}+n_{ξ}*n^{δ}*n_{λ}*n_{ε}*n_{ζ}*n^{λ}*n^{ξ}*n^{μ}*g^{νγ}+n_{ξ}*n^{γ}*n_{λ}*n^{ι}*n_{ι}*n^{λ}*n^{ξ}*n^{μ}*g^{νδ}*g_{εζ}+n_{ξ}*n^{δ}*n_{λ}*n^{ι}*n_{ι}*n^{λ}*n^{ξ}*n^{μ}*g^{νγ}*g_{εζ}+n_{λ}*n^{λ}*g^{νγ}*d_{ζ}^{δ}*d^{μ}_{ε}+n_{λ}*n^{λ}*d_{ζ}^{γ}*g^{νδ}*d^{μ}_{ε}+n_{λ}*n^{γ}*n^{ν}*n^{ο}*n^{λ}*n_{ο}*d_{ζ}^{δ}*d^{μ}_{ε}+n_{λ}*n^{δ}*n^{ν}*n^{ο}*n^{λ}*n_{ο}*d_{ζ}^{γ}*d^{μ}_{ε}+n_{λ}*n^{ο}*n^{λ}*n_{ο}*g^{νγ}*d_{ζ}^{δ}*d^{μ}_{ε}+n_{λ}*n^{ο}*n^{λ}*n_{ο}*d_{ζ}^{γ}*g^{νδ}*d^{μ}_{ε}+n_{ζ}*n^{γ}*n_{λ}*n^{ο}*n^{λ}*n_{ο}*g^{νδ}*d^{μ}_{ε}+n_{ζ}*n^{δ}*n_{λ}*n^{ο}*n^{λ}*n_{ο}*g^{νγ}*d^{μ}_{ε}+n_{λ}*n_{ξ}*n^{ν}*n^{ο}*n_{ζ}*n^{λ}*n^{ξ}*n_{ο}*g^{γδ}*d^{μ}_{ε}+n_{λ}*n_{ξ}*n^{γ}*n^{δ}*n^{ν}*n^{ο}*n_{ζ}*n^{λ}*n^{ξ}*n_{ο}*d^{μ}_{ε}+n_{ξ}*n^{γ}*n_{λ}*n^{ο}*n_{ζ}*n^{λ}*n^{ξ}*n_{ο}*g^{νδ}*d^{μ}_{ε}+n_{ξ}*n^{δ}*n_{λ}*n^{ο}*n_{ζ}*n^{λ}*n^{ξ}*n_{ο}*g^{νγ}*d^{μ}_{ε}+n_{λ}*n^{λ}*g^{νγ}*d_{ε}^{δ}*d^{μ}_{ζ}+n_{λ}*n^{λ}*d_{ε}^{γ}*g^{νδ}*d^{μ}_{ζ}+n_{λ}*n^{γ}*n^{ν}*n^{ο}*n^{λ}*n_{ο}*d_{ε}^{δ}*d^{μ}_{ζ}+n_{λ}*n^{δ}*n^{ν}*n^{ο}*n^{λ}*n_{ο}*d_{ε}^{γ}*d^{μ}_{ζ}+n_{λ}*n^{ο}*n^{λ}*n_{ο}*g^{νγ}*d_{ε}^{δ}*d^{μ}_{ζ}+n_{λ}*n^{ο}*n^{λ}*n_{ο}*d_{ε}^{γ}*g^{νδ}*d^{μ}_{ζ}+n_{ε}*n^{γ}*n_{λ}*n^{ο}*n^{λ}*n_{ο}*g^{νδ}*d^{μ}_{ζ}+n_{ε}*n^{δ}*n_{λ}*n^{ο}*n^{λ}*n_{ο}*g^{νγ}*d^{μ}_{ζ}+n_{λ}*n_{ξ}*n^{ν}*n^{ο}*n_{ε}*n^{λ}*n^{ξ}*n_{ο}*g^{γδ}*d^{μ}_{ζ}+n_{λ}*n_{ξ}*n^{γ}*n^{δ}*n^{ν}*n^{ο}*n_{ε}*n^{λ}*n^{ξ}*n_{ο}*d^{μ}_{ζ}+n_{ξ}*n^{γ}*n_{λ}*n^{ο}*n_{ε}*n^{λ}*n^{ξ}*n_{ο}*g^{νδ}*d^{μ}_{ζ}+n_{ξ}*n^{δ}*n_{λ}*n^{ο}*n_{ε}*n^{λ}*n^{ξ}*n_{ο}*g^{νγ}*d^{μ}_{ζ}");

    }
}
