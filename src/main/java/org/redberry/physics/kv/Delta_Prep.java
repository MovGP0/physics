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
package org.redberry.physics.kv;

import redberry.core.context.CC;
import redberry.core.tensor.Expression;
import redberry.core.tensor.SimpleTensor;
import redberry.core.tensor.Tensor;
import redberry.core.transformation.CalculateNumbers;
import redberry.core.transformation.IndexesInsertion;
import redberry.core.transformation.Transformation;
import redberry.core.transformation.Transformations;
import redberry.core.transformation.collect.CollecctEqualsInputPort;
import redberry.core.transformation.concurrent.EACScalars;
import redberry.core.transformation.concurrent.ExpandAndCollectTransformation;
import redberry.core.transformation.contractions.IndexesContractionsTransformation;
import redberry.core.transformation.substitutions.NaiveSubstitution;
import redberry.core.utils.Indicator;
import org.redberry.physics.util.SqrSubs;
import static org.redberry.physics.util.IndexesFactoryUtil.*;

/**
 *
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
public class Delta_Prep {
    public static void go(OneLoop loop) {
        evalD1(loop);
        evalD2(loop);
        evalD3(loop);
        evalD4(loop);
    }
    static final Transformation ec = new ExpandAndCollectTransformation(
            new CollecctEqualsInputPort(),
            Indicator.SYMBOL_INDICATOR,
            new Transformation[]{
                IndexesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                OneLoop.KRONECKER_DIMENSION.asSubstitution(),
                new SqrSubs((SimpleTensor) CC.parse("n_{\\alpha}")),
                CalculateNumbers.INSTANCE});

    static void evalD1(OneLoop loop) {
        for (Expression hK : loop.HATKs)
            loop.DELTA_1.eval(hK.asSubstitution());
        loop.DELTA_1.eval(ec, EACScalars.getTransformer(), CalculateNumbers.INSTANCE);
    }

    static void evalD2(OneLoop loop) {
        for (Expression hK : loop.HATKs)
            loop.DELTA_2.eval(hK.asSubstitution());
        loop.DELTA_2.eval(ec, EACScalars.getTransformer(), CalculateNumbers.INSTANCE);
    }

    static void evalD3(OneLoop loop) {
        Tensor[] hatkCombinations = CC.parse(
                "HATK^{\\mu\\nu}*HATK^{\\alpha}",
                "HATK^{\\alpha}*HATK^{\\mu\\nu}",
                "HATK^{\\mu}*HATK^{\\nu}*HATK^{\\alpha}");
        IndexesInsertion indexesInsertion = new IndexesInsertion(loop.matricesIndicator, createIndexes(hatkCombinations, "^{\\mu\\nu}_{\\alpha\\beta}"));

        for (int i = 0; i < hatkCombinations.length; ++i)
            hatkCombinations[i] = indexesInsertion.transform(hatkCombinations[i]);


        Expression[] hatkCombDone = new Expression[hatkCombinations.length];
        for (int i = 0; i < hatkCombinations.length; ++i) {

            Tensor hatkExpanded = hatkCombinations[i].clone();
            for (Expression hK : loop.HATKs)
                hK.asSubstitution().transform(hatkExpanded);
            hatkExpanded = MainScenario.smartEC(hatkExpanded, ec);
            hatkExpanded = Transformations.expandAndCollectAllScalars(hatkExpanded);
            hatkCombDone[i] = new Expression(hatkCombinations[i], hatkExpanded);
        }

        for (Expression t : hatkCombDone)
            NaiveSubstitution.create(t).transform(loop.DELTA_3);

        loop.DELTA_3.eval(ec, EACScalars.getTransformer(), CalculateNumbers.INSTANCE);

        int c1 = MainScenario.getElementsCount(loop.DELTA_3);
        for (Expression h : loop.HATKs)
            h.asSubstitution().transform(loop.DELTA_3);
        int c2 = MainScenario.getElementsCount(loop.DELTA_3);

        assert c1 == c2;
    }

    static void evalD4(OneLoop loop) {
        Tensor[] hatkCombinations = CC.parse(
                "HATK^{\\alpha\\beta}*HATK^{\\mu\\nu}",
                "HATK^{\\beta}*HATK^{\\alpha}*HATK^{\\mu}*HATK^{\\nu}",
                "HATK^{\\mu}*HATK^{\\nu}*HATK^{\\alpha\\beta}",
                "HATK^{\\mu}*HATK^{\\alpha\\beta}*HATK^{\\nu}",
                "HATK^{\\alpha\\beta}*HATK^{\\mu}*HATK^{\\nu}");
        IndexesInsertion indexesInsertion = new IndexesInsertion(loop.matricesIndicator, createIndexes(hatkCombinations, "^{\\mu\\nu}_{\\alpha\\beta}"));

        for (int i = 0; i < hatkCombinations.length; ++i)
            hatkCombinations[i] = indexesInsertion.transform(hatkCombinations[i]);

        Expression[] hatkCombDone = new Expression[hatkCombinations.length];
        for (int i = 0; i < hatkCombinations.length; ++i) {

            Tensor hatkExpanded = hatkCombinations[i].clone();
            for (Expression hK : loop.HATKs)
                hK.asSubstitution().transform(hatkExpanded);
            hatkExpanded = MainScenario.smartEC(hatkExpanded, ec);
            hatkExpanded = Transformations.expandAndCollectAllScalars(hatkExpanded);
            hatkCombDone[i] = new Expression(hatkCombinations[i], hatkExpanded);
        }

        for (Expression t : hatkCombDone)
            NaiveSubstitution.create(t).transform(loop.DELTA_4);

        loop.DELTA_4.eval(ec, EACScalars.getTransformer(), CalculateNumbers.INSTANCE);

        int c1 = MainScenario.getElementsCount(loop.DELTA_4);
        for (Expression h : loop.HATKs)
            h.asSubstitution().transform(loop.DELTA_4);
        int c2 = MainScenario.getElementsCount(loop.DELTA_4);

        assert c1 == c2;
    }
}
