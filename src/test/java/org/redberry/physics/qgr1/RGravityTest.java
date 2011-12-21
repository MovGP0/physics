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
package org.redberry.physics.qgr1;

import org.redberry.physics.qgr1.RGravity;
import org.junit.Test;
import redberry.core.tensor.Expression;
import redberry.core.transformation.CalculateNumbers;
import redberry.core.transformation.ExpandBrackets;
import redberry.core.transformation.Transformer;
import redberry.core.transformation.collect.CollectFactory;
import redberry.core.transformation.collect.CollectPowers;
import redberry.core.transformation.contractions.IndexesContractionsTransformation;

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
                IndexesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                rGravity.KRONECKER_DIMENSION.asSubstitution(),
                CollectFactory.createCollectAllEqualTerms(),
                CalculateNumbers.INSTANCE,
                new Transformer(CollectPowers.INSTANCE));
        System.out.println(loop);
    }
}
