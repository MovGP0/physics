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
package cc.redberry.physics.kv;

import cc.redberry.core.context.CC;
import cc.redberry.core.tensor.SimpleTensor;
import cc.redberry.core.tensor.Tensor;
import cc.redberry.core.utils.Indicator;
import cc.redberry.transformation.CalculateNumbers;
import cc.redberry.transformation.Transformation;
import cc.redberry.transformation.collect.CollecctEqualsInputPort;
import cc.redberry.transformation.concurrent.ExpandAndCollectTransformation;
import cc.redberry.transformation.contractions.IndicesContractionsTransformation;
import org.junit.Test;
import cc.redberry.physics.util.SqrSubs;



/**
 *
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
public class MainScenarioTest {
    public MainScenarioTest() {
    }

    @Test
    public void testSmartEC() {
        Tensor t = CC.parse("n_a*n_i*((a+b)*n_b+(c+d)*s_b)+c*(A_aib+A_aib)");
        Transformation ec = new ExpandAndCollectTransformation(
                new CollecctEqualsInputPort(),
                Indicator.SYMBOL_INDICATOR,
                new Transformation[]{
                    IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                    OneLoop.KRONECKER_DIMENSION.asSubstitution(),
                    new SqrSubs((SimpleTensor) CC.parse("n_{\\alpha}")),
                    CalculateNumbers.INSTANCE});

        t = MainScenario.smartEC(t, ec);
        System.out.println(t);
    }
    
    @Test
    public void testSmartEC1() {
        Tensor t = CC.parse("((-1/15)*LAMBDA*LAMBDA+1/45*LAMBDA*LAMBDA+2/45*LAMBDA*LAMBDA)*(1/2*n_{b}*d_{r}^{e}*d_{s}^{k}+1/2*n_{b}*d_{r}^{k}*d_{s}^{e})*(G^brs_ek+G^bsr_ek)");
        Transformation ec = new ExpandAndCollectTransformation(
                new CollecctEqualsInputPort(),
                Indicator.SYMBOL_INDICATOR,
                new Transformation[]{
                    IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                    OneLoop.KRONECKER_DIMENSION.asSubstitution(),
                    new SqrSubs((SimpleTensor) CC.parse("n_{\\alpha}")),
                    CalculateNumbers.INSTANCE});

        t = MainScenario.smartEC(t, ec);
        System.out.println(t);
    }
}
