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
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package redberryphysics.core.oneloop;

import org.junit.Ignore;
import redberry.core.tensor.Product;
import redberry.core.tensor.Sum;
import redberry.core.transformation.Transformation;
import redberry.core.transformation.collect.CollectTerms;
import redberry.core.transformation.collect.EqualsSplitCriteria;
import redberry.core.transformation.substitutions.TensorTreeIndicatorImpl;
import redberry.core.tensor.Tensor;
import java.util.Deque;
import redberry.core.transformation.collect.CollectPatternParser;
import org.junit.Test;
import redberry.core.context.CC;
import redberry.core.context.ToStringMode;
import redberry.core.parser.ParserIndexes;
import redberry.core.tensor.MultiTensor;
import redberry.core.tensor.TensorIterator;
import redberry.core.tensor.test.TTest;
import redberry.core.transformation.IndexesInsertion;
import redberry.core.transformation.Transformations;
import redberry.core.transformation.collect.PatternCollectManager;
import redberry.core.transformation.collect.SplitPattern;
import redberry.core.transformation.substitutions.SubstitutionsFactory;
import redberryphysics.core.oneloop.OneLoop;
import static redberryphysics.core.oneloop.OneLoop.*;
import static core.TAssert.*;
import static org.junit.Assert.*;

/**
 *
 * @author stas
 */
public class OneLoopTest {
    @Ignore
    @Test
    public void substitutions1() {
        OneLoop loop = new OneLoop();
        RR = loop.RR_SUBSTITUTION.transform(RR);
        RR = loop.DELTA_1_SUBSTITUTION.transform(RR);
        RR = loop.DELTA_2_SUBSTITUTION.transform(RR);
        RR = loop.DELTA_3_SUBSTITUTION.transform(RR);
        RR = loop.DELTA_4_SUBSTITUTION.transform(RR);
        RR = Transformations.multiplyNumbers(RR);
        RR = Transformations.sumNumbers(RR);
        Tensor copy = RR.clone();
        assertTrue(TTest.testIsScalar(RR));
        assertIndexes(RR);
        System.out.println(((MultiTensor) RR).size());
        RR = Transformations.expandBrackets(RR);
        assertTrue(TTest.testIsScalar(RR));
        assertIndexes(RR);
        System.out.println(((MultiTensor) RR).size());
        Deque<SplitPattern> deque = CollectPatternParser.parse("(HATK_{\\mu}"
                + "|HATK_{\\mu\\nu}"
                + "|HATK_{\\mu\\nu\\alpha}"
                + "|HATK_{\\mu\\nu\\alpha\\beta}"
                + "|n_{\\mu}"
                + "|(R_{\\mu\\nu}|R_{\\mu\\nu\\alpha\\beta}))");
        PatternCollectManager collectManager = new PatternCollectManager(deque);

        RR = collectManager.transform(RR);
        assertTrue(TTest.testIsScalar(RR));
        assertIndexes(RR);
        System.out.println(((MultiTensor) RR).size());
        RR = Transformations.multiplyNumbers(RR);
        RR = Transformations.sumNumbers(RR);
    }
 @Ignore
    @Test
    public void substitutions2_L_TO_2() {
        OneLoop loop = new OneLoop();
        RR = loop.RR_SUBSTITUTION.transform(RR);
        RR = loop.DELTA_1_SUBSTITUTION.transform(RR);
        RR = loop.DELTA_2_SUBSTITUTION.transform(RR);
        RR = loop.DELTA_3_SUBSTITUTION.transform(RR);
        RR = loop.DELTA_4_SUBSTITUTION.transform(RR);
        RR = SubstitutionsFactory.createSubstitution("L=2").transform(RR);
        RR = Transformations.calculateNumbers(RR);
        assertTrue(TTest.testIsScalar(RR));
        assertIndexes(RR);
        System.out.println(RR);
    }

    public static void cutOfOddN(Sum sum) {
        TensorIterator iterator = sum.iterator();
        Tensor p;
        while (iterator.hasNext()) {
            p = iterator.next();
            if (!(p instanceof Product))
                continue;
            int count = 0;
            for (Tensor t : p)
                if (TTest.testEqualstensorStructure(t, OneLoop.N))
                    count++;
            if (count % 2 != 0)
                iterator.remove();
        }
    }

    @Test
    public void substitutions3_MATRICES() {
        OneLoop loop = new OneLoop();

        RR = loop.RR_SUBSTITUTION.transform(RR);
        RR = loop.DELTA_1_SUBSTITUTION.transform(RR);
        RR = loop.DELTA_2_SUBSTITUTION.transform(RR);
        RR = loop.DELTA_3_SUBSTITUTION.transform(RR);
        RR = loop.DELTA_4_SUBSTITUTION.transform(RR);
        RR = SubstitutionsFactory.createSubstitution("L=2").transform(RR);
        RR = Transformations.calculateNumbers(RR);

        IndexesInsertion indexesInsertion = new IndexesInsertion(new TensorTreeIndicatorImpl(matrices), ParserIndexes.parse("^ab_ab"));
        RR = indexesInsertion.transform(RR);
        assertIndexes(RR);

        RR = loop.MATRIX_HATK_1_SUBSTITUTION.transform(RR);
        RR = loop.MATRIX_HATK_2_SUBSTITUTION.transform(RR);
//      
//        RR = Transformations.renameConflictingIndexes(RR);
        RR = Transformations.expandBrackets(RR);
        System.out.println(((MultiTensor) RR).size());
        RR = Transformations.calculateNumbers(RR);
        RR = Transformations.contractMetrics(RR);
        RR = SubstitutionsFactory.createSubstitution("d^a_a=4").transform(RR);
        RR = Transformations.calculateNumbers(RR);
        assertIndexes(RR);
        Transformation collect = new CollectTerms(EqualsSplitCriteria.INSTANCE);
        RR = collect.transform(RR);
        System.out.println(((MultiTensor) RR).size());
        System.out.println(RR.toString(ToStringMode.UTF8));
    }
    
    static Tensor TEMP = CC.parse("(-48/5)*n^{\\delta }*n^{\\gamma }*n_{\\sigma }*n_{\\lambda }*R^{\\lambda \\nu }_{\\gamma \\nu }*R^{\\sigma \\beta }_{\\delta \\beta }+(-21/5)*n^{\\delta }*n^{\\beta }*n_{\\sigma }*n_{\\lambda }*R^{\\lambda \\nu \\alpha }_{\\nu }*R^{\\sigma }_{\\alpha \\delta \\beta }+(-2/5)*n^{\\delta }*n^{\\alpha }*n_{\\sigma }*n_{\\lambda }*R^{\\lambda \\nu \\beta }_{\\nu }*R^{\\sigma }_{\\alpha \\delta \\beta }+(-63/5)*n^{\\delta }*n^{\\alpha }*n^{\\beta }*n^{\\gamma }*n_{\\sigma }*n_{\\lambda }*R^{\\lambda \\nu }_{\\gamma \\nu }*R^{\\sigma }_{\\alpha \\delta \\beta }+8/5*n^{\\delta }*n^{\\gamma }*n_{\\sigma }*n_{\\lambda }*R^{\\lambda }_{\\delta }^{\\beta \\mu }*R^{\\sigma }_{\\beta \\mu \\gamma }+8/5*n^{\\delta }*n^{\\beta }*n_{\\sigma }*n_{\\lambda }*R^{\\lambda }_{\\delta }^{\\gamma \\mu }*R^{\\sigma }_{\\beta \\mu \\gamma }+6*n^{\\delta }*n^{\\alpha }*n^{\\beta }*n^{\\gamma }*n_{\\sigma }*n_{\\lambda }*R^{\\lambda }_{\\delta \\alpha }^{\\mu }*R^{\\sigma }_{\\beta \\mu \\gamma }+2/5*n^{\\delta }*n^{\\gamma }*n_{\\sigma }*n_{\\lambda }*R^{\\lambda \\alpha }_{\\delta }^{\\mu }*R^{\\sigma }_{\\gamma \\mu \\alpha }+2/5*n^{\\delta }*n^{\\alpha }*n_{\\sigma }*n_{\\lambda }*R^{\\lambda \\gamma }_{\\delta }^{\\mu }*R^{\\sigma }_{\\gamma \\mu \\alpha }+5/2*n_{\\sigma }*n_{\\lambda }*R^{\\lambda \\nu \\alpha }_{\\nu }*R^{\\sigma \\gamma }_{\\alpha \\gamma }+(-1/5)*n_{\\sigma }*n_{\\lambda }*R^{\\lambda \\beta \\mu }_{\\beta }*R^{\\sigma \\delta }_{\\delta \\mu }+4/5*n^{\\alpha }*n^{\\beta }*n_{\\sigma }*n_{\\lambda }*R^{\\lambda }_{\\alpha }^{\\mu }_{\\beta }*R^{\\sigma \\delta }_{\\delta \\mu }+(-23/5)*n^{\\alpha }*n^{\\beta }*n_{\\sigma }*n_{\\lambda }*R^{\\lambda \\nu \\delta }_{\\nu }*R^{\\sigma }_{\\alpha \\delta \\beta }+9/20*n_{\\sigma }*n_{\\lambda }*R^{\\lambda \\gamma \\alpha \\mu }*R^{\\sigma }_{\\gamma \\alpha \\mu }+(-3/5)*n^{\\alpha }*n^{\\beta }*n_{\\sigma }*n_{\\lambda }*R^{\\lambda \\gamma }_{\\beta }^{\\mu }*R^{\\sigma }_{\\gamma \\alpha \\mu }+(-1)*n^{\\alpha }*n^{\\beta }*n_{\\sigma }*n_{\\lambda }*R^{\\lambda }_{\\beta }^{\\delta \\mu }*R^{\\sigma }_{\\alpha \\delta \\mu }+6/5*n_{\\sigma }*n_{\\lambda }*R^{\\lambda \\gamma \\alpha \\mu }*R^{\\sigma }_{\\alpha \\gamma \\mu }+(-23/5)*n^{\\alpha }*n^{\\beta }*n_{\\sigma }*n_{\\lambda }*R^{\\lambda \\gamma }_{\\beta }^{\\mu }*R^{\\sigma }_{\\alpha \\gamma \\mu }+7/5*n^{\\gamma }*n_{\\sigma }*n_{\\lambda }*R^{\\lambda \\alpha }_{\\gamma }^{\\mu }*R^{\\sigma }_{\\mu \\alpha \\delta }*n^{\\delta }+7/5*n^{\\beta }*n_{\\sigma }*n_{\\lambda }*R^{\\lambda }_{\\beta }^{\\alpha \\mu }*R^{\\sigma }_{\\mu \\alpha \\delta }*n^{\\delta }+7/5*n^{\\alpha }*n_{\\sigma }*n_{\\lambda }*R^{\\lambda \\gamma }_{\\gamma }^{\\mu }*R^{\\sigma }_{\\mu \\alpha \\delta }*n^{\\delta }+21/5*n^{\\alpha }*n^{\\beta }*n^{\\gamma }*n_{\\sigma }*n_{\\lambda }*R^{\\lambda }_{\\beta \\gamma }^{\\mu }*R^{\\sigma }_{\\mu \\alpha \\delta }*n^{\\delta }+(-9/5)*n^{\\alpha }*n_{\\sigma }*n_{\\lambda }*R^{\\lambda \\gamma }_{\\gamma \\delta }*R^{\\sigma \\nu }_{\\alpha \\nu }*n^{\\delta }+1/5*n^{\\gamma }*n_{\\sigma }*n_{\\lambda }*R^{\\lambda }_{\\delta \\gamma }^{\\mu }*R^{\\sigma \\beta }_{\\beta \\mu }*n^{\\delta }+3/5*n^{\\alpha }*n^{\\beta }*n^{\\gamma }*n_{\\sigma }*n_{\\lambda }*R^{\\lambda }_{\\delta \\gamma }^{\\mu }*R^{\\sigma }_{\\alpha \\beta \\mu }*n^{\\delta }+12/5*n^{\\mu }*n^{\\gamma }*{\\nu}_{\\lambda }*R^{\\alpha }_{\\gamma }*R^{\\lambda }_{\\nu \\alpha \\mu }*n^{\\nu }+3/5*n^{\\mu }*n^{\\beta }*{\\nu}_{\\lambda }*R_{\\beta }^{\\alpha }*R^{\\lambda }_{\\nu \\alpha \\mu }*n^{\\nu }+3/5*n^{\\mu }*n^{\\alpha }*{\\nu}_{\\lambda }*R^{\\gamma }_{\\gamma }*R^{\\lambda }_{\\nu \\alpha \\mu }*n^{\\nu }+36/5*n^{\\mu }*n^{\\alpha }*n^{\\beta }*n^{\\gamma }*{\\nu}_{\\lambda }*R_{\\beta \\gamma }*R^{\\lambda }_{\\nu \\alpha \\mu }*n^{\\nu }+9/5*n^{\\mu }*n^{\\gamma }*{\\nu}_{\\lambda }*R^{\\beta }_{\\mu }*R^{\\lambda }_{\\beta \\gamma \\nu }*n^{\\nu }+9/5*n^{\\mu }*n^{\\alpha }*{\\nu}_{\\lambda }*R_{\\alpha \\mu }*R^{\\lambda \\gamma }_{\\gamma \\nu }*n^{\\nu }+3/5*n^{\\mu }*n^{\\gamma }*{\\nu}_{\\lambda }*R^{\\sigma \\alpha }_{\\gamma \\mu }*R^{\\lambda }_{\\nu \\alpha \\sigma }*n^{\\nu }+3/5*n^{\\mu }*n^{\\beta }*{\\nu}_{\\lambda }*R^{\\sigma }_{\\beta }^{\\alpha }_{\\mu }*R^{\\lambda }_{\\nu \\alpha \\sigma }*n^{\\nu }+(-3/5)*n^{\\mu }*n^{\\alpha }*{\\nu}_{\\lambda }*R^{\\sigma \\gamma }_{\\gamma \\mu }*R^{\\lambda }_{\\nu \\alpha \\sigma }*n^{\\nu }+(-9/5)*n^{\\mu }*n^{\\alpha }*n^{\\beta }*n^{\\gamma }*{\\nu}_{\\lambda }*R^{\\sigma }_{\\beta \\gamma \\mu }*R^{\\lambda }_{\\nu \\alpha \\sigma }*n^{\\nu }+(-6/5)*n^{\\mu }*n^{\\beta }*{\\nu}_{\\lambda }*R^{\\sigma \\gamma }_{\\beta \\mu }*R^{\\lambda }_{\\gamma \\nu \\sigma }*n^{\\nu }+(-6/5)*n^{\\mu }*n^{\\alpha }*{\\nu}_{\\lambda }*R^{\\sigma }_{\\alpha }^{\\gamma }_{\\mu }*R^{\\lambda }_{\\gamma \\nu \\sigma }*n^{\\nu }+(-9/5)*n^{\\gamma }*n_{\\lambda }*R^{\\beta \\mu }*R^{\\lambda }_{\\gamma \\beta \\mu }+36/5*n^{\\gamma }*n^{\\alpha }*n^{\\beta }*n_{\\lambda }*R_{\\alpha }^{\\mu }*R^{\\lambda }_{\\gamma \\beta \\mu }+(-36/5)*n^{\\gamma }*n_{\\lambda }*R^{\\beta }_{\\gamma }*R^{\\lambda \\nu }_{\\beta \\nu }+108/5*n^{\\gamma }*n^{\\alpha }*n^{\\beta }*n_{\\lambda }*R_{\\alpha \\gamma }*R^{\\lambda \\nu }_{\\beta \\nu }+(-57/10)*n^{\\gamma }*n_{\\lambda }*R^{\\beta }_{\\beta }*R^{\\lambda \\nu }_{\\gamma \\nu }+(-9/5)*n^{\\gamma }*n_{\\lambda }*R^{\\sigma \\beta \\mu }_{\\gamma }*R^{\\lambda }_{\\sigma \\beta \\mu }+36/5*n^{\\gamma }*n^{\\alpha }*n^{\\beta }*n_{\\lambda }*R^{\\sigma }_{\\alpha }^{\\mu }_{\\gamma }*R^{\\lambda }_{\\sigma \\beta \\mu }+3/5*n^{\\gamma }*n_{\\lambda }*R^{\\sigma \\nu \\beta }_{\\nu }*R^{\\lambda }_{\\beta \\sigma \\gamma }+(-12/5)*n^{\\gamma }*n^{\\alpha }*n^{\\beta }*n_{\\lambda }*R^{\\sigma \\nu }_{\\alpha \\nu }*R^{\\lambda }_{\\beta \\sigma \\gamma }+(-18/5)*n^{\\gamma }*n_{\\lambda }*R^{\\sigma \\beta }_{\\beta \\gamma }*R^{\\lambda \\nu }_{\\sigma \\nu }+72/5*n^{\\gamma }*n^{\\alpha }*n^{\\beta }*n_{\\lambda }*R^{\\sigma }_{\\alpha \\beta \\gamma }*R^{\\lambda \\nu }_{\\sigma \\nu }+3*n^{\\gamma }*n_{\\lambda }*R^{\\sigma \\beta }_{\\beta }^{\\mu }*R^{\\lambda }_{\\mu \\sigma \\gamma }+(-12)*n^{\\gamma }*n^{\\alpha }*n^{\\beta }*n_{\\lambda }*R^{\\sigma }_{\\alpha \\beta }^{\\mu }*R^{\\lambda }_{\\mu \\sigma \\gamma }+(-6/5)*n^{\\gamma }*n_{\\lambda }*R^{\\beta \\nu }*R^{\\lambda }_{\\beta \\nu \\gamma }+18/5*n^{\\gamma }*n_{\\lambda }*R^{\\beta \\nu }*R^{\\lambda }_{\\gamma \\nu \\beta }+(-3/5)*n^{\\gamma }*n_{\\lambda }*R^{\\beta }_{\\beta }*R^{\\lambda \\nu }_{\\nu \\gamma }+(-18/5)*n^{\\gamma }*n_{\\lambda }*R^{\\sigma }_{\\gamma }^{\\nu \\beta }*R^{\\lambda }_{\\nu \\sigma \\beta }+(-21/5)*n^{\\gamma }*n_{\\lambda }*R^{\\sigma \\beta \\nu }_{\\gamma }*R^{\\lambda }_{\\beta \\nu \\sigma }+9/5*n^{\\gamma }*n_{\\lambda }*R^{\\sigma \\beta \\nu }_{\\gamma }*R^{\\lambda }_{\\sigma \\nu \\beta }+(-18/5)*n^{\\gamma }*n_{\\lambda }*R^{\\sigma \\beta \\nu }_{\\beta }*R^{\\lambda }_{\\sigma \\nu \\gamma }+39/5*n^{\\gamma }*n_{\\lambda }*R^{\\sigma \\beta \\nu }_{\\beta }*R^{\\lambda }_{\\gamma \\nu \\sigma }+9/5*n^{\\gamma }*n_{\\lambda }*R^{\\sigma \\nu \\alpha }_{\\nu }*R^{\\lambda }_{\\gamma \\sigma \\alpha }+(-9/5)*n_{\\lambda }*R^{\\sigma \\nu }_{\\nu }^{\\alpha }*R^{\\lambda }_{\\gamma \\sigma \\alpha }*n^{\\gamma }+36/5*n^{\\mu }*n^{\\nu }*n_{\\lambda }*R^{\\sigma }_{\\mu \\nu }^{\\alpha }*R^{\\lambda }_{\\gamma \\sigma \\alpha }*n^{\\gamma }+(-84/5)*n^{\\mu }*n^{\\nu }*n_{\\lambda }*R^{\\sigma \\alpha }_{\\mu \\alpha }*R^{\\lambda }_{\\gamma \\nu \\sigma }*n^{\\gamma }+36/5*n^{\\mu }*n^{\\nu }*n_{\\lambda }*R^{\\sigma \\alpha }_{\\mu \\alpha }*R^{\\lambda }_{\\sigma \\nu \\gamma }*n^{\\gamma }+72/5*n^{\\mu }*n^{\\nu }*n_{\\lambda }*R^{\\sigma }_{\\mu }^{\\alpha }_{\\gamma }*R^{\\lambda }_{\\nu \\alpha \\sigma }*n^{\\gamma }+(-3/5)*n_{\\lambda }*R^{\\sigma \\nu \\alpha }_{\\gamma }*R^{\\lambda }_{\\alpha \\nu \\sigma }*n^{\\gamma }+12/5*n^{\\mu }*n^{\\nu }*n_{\\lambda }*R^{\\sigma }_{\\mu }^{\\alpha }_{\\gamma }*R^{\\lambda }_{\\alpha \\nu \\sigma }*n^{\\gamma }+42/5*n^{\\mu }*n^{\\nu }*n_{\\lambda }*R^{\\beta }_{\\beta }*R^{\\lambda }_{\\nu \\gamma \\mu }*n^{\\gamma }+(-33/5)*n_{\\lambda }*R^{\\alpha \\nu }*R^{\\lambda }_{\\nu \\alpha \\gamma }*n^{\\gamma }+132/5*n^{\\mu }*n^{\\nu }*n_{\\lambda }*R^{\\alpha }_{\\mu }*R^{\\lambda }_{\\nu \\alpha \\gamma }*n^{\\gamma }");
        
    @Test
    public void substitutions() {
        System.out.println(((Sum)TEMP).size());
        cutOfOddN((Sum)TEMP);
        System.out.println(((Sum)TEMP).size());        
    }
}
