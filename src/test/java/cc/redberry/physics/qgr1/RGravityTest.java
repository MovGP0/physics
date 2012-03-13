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
package cc.redberry.physics.qgr1;

import cc.redberry.core.tensor.Expression;
import cc.redberry.transformation.CalculateNumbers;
import cc.redberry.transformation.ExpandBrackets;
import cc.redberry.transformation.Transformer;
import cc.redberry.transformation.collect.CollectFactory;
import cc.redberry.transformation.collect.CollectPowers;
import cc.redberry.transformation.contractions.IndicesContractionsTransformation;
import org.junit.Test;



/**
 *
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
public class RGravityTest {
    public RGravityTest() {
    }

    @Test
    public void test2G2Gr() {
        RGravity rGravity = new RGravity();
        System.out.println(rGravity.SCALAR2_GRAVITON2_VERTEX);
    }

    @Test
    public void oneLoop1() {
        RGravity rGravity = new RGravity();
        Expression loop = new Expression("G1^mnab[k_a] = T^mn_rs[p_a,p_a+k_a]*G[p_a+k_a]*G_[p_a]*T^rsab[p_a,p_a+k_a]");
        rGravity.substituteAllToExpression(loop);
        loop.eval(
                new Transformer(ExpandBrackets.EXPAND_EXCEPT_SYMBOLS),
                IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                rGravity.KRONECKER_DIMENSION.asSubstitution(),
                CollectFactory.createCollectAllEqualTerms(),
                CalculateNumbers.INSTANCE,
                new Transformer(CollectPowers.INSTANCE));
        System.out.println(loop);
    }
}
