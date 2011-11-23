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

import redberry.core.transformation.concurrent.EACScalars;
import redberry.core.context.CC;
import redberry.core.tensor.Tensor;
import org.junit.Test;
import redberry.core.context.ToStringMode;
import redberry.core.tensor.Expression;
import redberry.core.tensor.SimpleTensor;
import redberry.core.tensor.Sum;
import redberry.core.transformation.CalculateNumbers;
import redberry.core.transformation.Transformation;
import redberry.core.transformation.Transformations;
import redberry.core.transformation.Transformer;
import redberry.core.transformation.collect.CollecctEqualsInputPort;
import redberry.core.transformation.collect.CollectFactory;
import redberry.core.transformation.collect.CollectPowers;
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
public class Delta_PrepTest {
    public Delta_PrepTest() {
    }

    @Test
    public void test() {
        Tensor t = CC.parse("(-2)*(3/8*b+1/4*b*beta)*n_{e}*n^{g}*n^{m}*d_{c}^{d}+(-2)*(3/16*b+3/4+1/8*b*beta+1/2*beta)*n_{e}*n^{d}*n^{m}*d_{c}^{g}");
        t = Delta_Prep.ec.transform(t);
        t = EACScalars.getTransformer().transform(t);
        t = CalculateNumbers.INSTANCE.transform(t);
        System.out.println(t);
    }

    @Test
    public void testEvalD1() {
        OneLoop loop = new OneLoop();
        loop.insertIndexes();
        loop.substituteL();
        loop.evalHatK();
        Delta_Prep.evalD1(loop);
        loop.DELTA_1.asSubstitution();
        System.out.println(loop.DELTA_1.toString(ToStringMode.UTF8));
    }

    @Test
    public void testEvalD2() {
        OneLoop loop = new OneLoop();
        loop.insertIndexes();
        loop.substituteL();
        loop.evalHatK();
        Delta_Prep.evalD2(loop);
        loop.DELTA_2.asSubstitution();
        System.out.println(loop.DELTA_2.toString(ToStringMode.UTF8));
    }

    @Test
    public void testEvalD3() {
        OneLoop loop = new OneLoop();
        loop.insertIndexes();
        loop.substituteL();
        loop.evalHatK();
        Delta_Prep.evalD3(loop);
        loop.DELTA_3.asSubstitution();
        System.out.println(loop.DELTA_3.toString(ToStringMode.UTF8));
    }

    @Test
    public void testEvalD4() {
        OneLoop loop = new OneLoop();
        loop.insertIndexes();
        loop.substituteL();
        loop.evalHatK();
        Delta_Prep.evalD4(loop);
        loop.DELTA_4.asSubstitution();
        System.out.println(loop.DELTA_4.toString(ToStringMode.UTF8));
    }

    @Test
    public void testEvalD3RR() {
        OneLoop loop = new OneLoop();
        loop.insertIndexes();
        loop.substituteL();
        loop.evalHatK();
        Delta_Prep.evalD3(loop);

        System.out.println(loop.DELTA_3.toString(ToStringMode.UTF8));
        Transformation ec = new ExpandAndCollectTransformation(
                new CollecctEqualsInputPort(),
                Indicator.SYMBOL_INDICATOR,
                new Transformation[]{
                    IndexesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                    loop.KRONECKER_DIMENSION.asSubstitution(),
                    new SqrSubs((SimpleTensor) CC.parse("n_{\\alpha}")),
                    CalculateNumbers.INSTANCE});
        long startTime = System.currentTimeMillis();
        System.out.println("Sustitution to FIRST");


        Tensor firstSummand = ((Sum) loop.RR.right()).getElements().get(1);

        //System.out.println(firstSummand.toString(ToStringMode.UTF8));

        firstSummand = loop.L.asSubstitution().transform(firstSummand);
        firstSummand = CalculateNumbers.INSTANCE.transform(firstSummand);
        firstSummand = loop.RIMAN.asSubstitution().transform(firstSummand);
        firstSummand = loop.RICCI.asSubstitution().transform(firstSummand);
        firstSummand = Transformations.expandBrackets(firstSummand);
        firstSummand = Transformations.contractMetrics(firstSummand);
        firstSummand = CollectFactory.createCollectEqualTerms1().transform(firstSummand);
        firstSummand = CalculateNumbers.INSTANCE.transform(firstSummand);

        for (Expression h : loop.HATKs)
            h.asSubstitution().transform(firstSummand);

        loop.DELTA_3.asSubstitution().transform(firstSummand);

        firstSummand = Transformations.contractMetrics(firstSummand);
        firstSummand = loop.KRONECKER_DIMENSION.asSubstitution().transform(firstSummand);
        firstSummand = CalculateNumbers.INSTANCE.transform(firstSummand);

        System.out.println(((Sum) loop.DELTA_3.right()).size());
        System.out.println(((Sum) loop.HATK_1.right()).size());

        System.out.println("First elements: " + MainScenario.getElementsCount(firstSummand));

        System.out.println("Exp+Collect");

        firstSummand = MainScenario.smartEC(firstSummand, ec);

        System.out.println(firstSummand);
        System.out.println("Collect Scalar");
        Transformation sc = new ExpandAndCollectTransformation(ScalarsSplitCriteria.INSTANCE,
                Indicator.FALSE_INDICATOR, new Transformation[]{CalculateNumbers.INSTANCE});
        firstSummand = sc.transform(firstSummand);
        firstSummand = CalculateNumbers.INSTANCE.transform(firstSummand);
        long stopTime = System.currentTimeMillis();


        System.out.println(firstSummand);
        System.out.println("Done. Elements: " + MainScenario.getElementsCount(firstSummand));

        System.out.println("First term time = " + (stopTime - startTime) + "ms");
        System.out.println("Total RR tertms count " + ((Sum) loop.RR.right()).size());
        firstSummand = new Transformer(CollectPowers.INSTANCE).transform(firstSummand);
        System.out.println(firstSummand);

    }
}
