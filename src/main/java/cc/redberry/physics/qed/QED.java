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
package cc.redberry.physics.qed;

import cc.redberry.core.context.CC;
import cc.redberry.core.tensor.Expression;
import cc.redberry.transformation.CalculateNumbers;
import cc.redberry.transformation.ExpandBrackets;
import cc.redberry.transformation.Transformer;
import cc.redberry.transformation.contractions.IndicesContractionsTransformation;
import cc.redberry.transformation.integral.CollectIntegralFromSum;
import cc.redberry.transformation.integral.ExpandIntegral;
import cc.redberry.physics.ToFourier;



/**
 *
 * @author Stanislav Poslavsky
 */
public class QED {
    public final Expression Action =
            new Expression("S = Integral[Lagrangian[x_a],x_a]");
    public final Expression Lagrangian =
            new Expression("Lagrangian[x_a] = (-1/4)*F_{mn}[x_a]*F^{mn}[x_a]");
    public final Expression F =
            new Expression("F_{mn}[x_a] = D[A_{m}[x_a],x^n]-D[A_{n}[x_a],x^m]");
    final Expression[] initial = new Expression[]{Action, Lagrangian, F};

    public QED() {
        Action.eval(
                Lagrangian,
                F,
                new Transformer(ExpandBrackets.EXPAND_ALL),
                new Transformer(ExpandIntegral.INSTANCE),
                new Transformer(new ToFourier(CC.parseSimple("x_a"), CC.parseSimple("p_a"))),
                IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                CalculateNumbers.INSTANCE,
                CollectIntegralFromSum.INSTANCE);
    }
}
