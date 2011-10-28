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

import redberry.core.tensor.Tensor;
import redberry.core.tensor.Sum;
import org.junit.BeforeClass;
import org.junit.Test;
import redberry.core.context.ToStringMode;
import redberry.core.tensor.MultiTensor;
import redberry.core.transformation.CalculateNumbers;
import redberry.core.transformation.ExpandBrackets;
import redberry.core.transformation.IndexesInsertion;
import redberry.core.transformation.Transformer;
import redberry.core.transformation.collect.CollectFactory;
import redberry.core.transformation.contractions.IndexesContractionsTransformation;
import redberryphysics.core.util.IndexesFactoryUtil;
import redberryphysics.core.util.ParallelCollect;
import static core.TAssert.*;

/**
 *
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 * @author Konstantin Kiselev
 */
public class OneLoop1Test {
    public OneLoop1Test() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
    }
//
//    @Test
//    public void test() {
//        OneLoop1 loop1 = new OneLoop1();
//        System.out.println(loop1.MATRIX_K.toString(ToStringMode.UTF8));
//        System.out.println(loop1.MATRIX_K_INV.toString(ToStringMode.UTF8));
//    }
//
//    @Test
//    public void test1() {
//        OneLoop1 loop1 = new OneLoop1();
//        System.out.println(loop1.HATK_1.toString(ToStringMode.UTF8));
//        System.out.println(loop1.HATK_2.toString(ToStringMode.UTF8));
//        System.out.println(loop1.DELTA_1.toString(ToStringMode.UTF8));
//        System.out.println(loop1.DELTA_2.toString(ToStringMode.UTF8));
//    }

    @Test
    public void test2() {
        OneLoop1 loop1 = new OneLoop1();
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
                new Transformer(ExpandBrackets.INSTANCE));
        System.out.println("Expand brackets ... done");

        System.out.println("Indexes contractions ...");
        loop1.DELTA_2.eval(
                IndexesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                loop1.KRONECKER_DIMENSION.asSubstitution(),
                CalculateNumbers.INSTANCE);
        System.out.println("Indexes contractions ... done");
        assertIndexes(loop1.DELTA_2.right());
        System.out.println("Collecting " + ((MultiTensor) loop1.DELTA_2.right()).size() + " elements ...");
        Sum sum = new Sum();
        int count = 0;
        for (Tensor t : loop1.DELTA_2.right()) {
            sum.add(t);
            if (count++ >= 1000)
                break;
        }
        loop1.DELTA_2.setRight(sum);
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
    }
}
