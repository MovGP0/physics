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
public class RGravity {
    public final Expression KRONECKER_DIMENSION =
            new Expression("d^{m}_{m} = 4");
    public final Expression P =
            new Expression("P^{mnab} = (1/2)*(g^{ma}*g^{nb}+g^{mb}*g^{na}-g^{mn}*g^{ab}) ");
    public final Expression I =
            new Expression("I^mnab = (1/2)*(g^ma*g^nb+g^mb*g^na)");
    /*
     * fields which can be evalueted
     */
    public final Expression GRAVITON_PROPAGATOR =
            new Expression("G^{mnab}[k_n] = (I*P^{mnab})/(k_n*k^n)");
    public final Expression SCALAR_PROPAGATOR =
            new Expression("G[k_n] = I/(k_n*k^n)");
    public final Expression SCALAR2_GRAVITON_VERTEX =
            new Expression("T^{mn}[k^n,p^n] = -(I/2)*G*(k^n*p^m+k^m*p^n-g^{mn}*(k_a*p^a-m*m))");
    public final Expression SCALAR2_GRAVITON2_VERTEX =
            new Expression("T^{nlrs}[p^n,k^n] = I*G*G*( (I^nlad*I^rsb_d-(1/4)*(g^nl*I^rsab+g^rs*I^nlab))*(p_a*k_b+p_b*k_a)-"
            + "(1/2)*(I^nlrs-(1/2)*g^nl*g^rs)*(p_a*k^a-m*m))");
    public final Expression[] rules = new Expression[]{
        GRAVITON_PROPAGATOR,
        SCALAR_PROPAGATOR,
        SCALAR2_GRAVITON2_VERTEX,
        SCALAR2_GRAVITON_VERTEX};

    public RGravity() {
        evalG_propagator();
        evalS2G2_vertex();
    }

    private void evalG_propagator() {
        GRAVITON_PROPAGATOR.eval(P.asSubstitution());
    }

    private void evalS2G2_vertex() {
        SCALAR2_GRAVITON2_VERTEX.eval(
                I.asSubstitution(),
                new Transformer(ExpandBrackets.EXPAND_EXCEPT_SYMBOLS),
                IndexesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                KRONECKER_DIMENSION.asSubstitution(),
                CollectFactory.createCollectAllEqualTerms(),
                CalculateNumbers.INSTANCE,
                new Transformer(CollectPowers.INSTANCE));
    }

    public void substituteAllToExpression(Expression e) {
        for (Expression r : rules)
            e.eval(r.asSubstitution());
    }
}
