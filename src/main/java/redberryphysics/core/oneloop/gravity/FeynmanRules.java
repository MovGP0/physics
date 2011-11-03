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
package redberryphysics.core.oneloop.gravity;

import redberry.core.tensor.Expression;

/**
 *
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
public class FeynmanRules {
    public final Expression KRONECKER_DIMENSION =
            new Expression("d^{m}_{m} = 4");
    public final Expression P =
            new Expression("P^{mnab} = (1/2)*(g^{ma}*g^{nb}+g^{mb}*g^{na}-g^{mn}*g^{ab}) ");
    public final Expression propagator =
            new Expression("G^{mnab}[k_n] = (I*P^{mnab})/(k_n*k^n)");
    public final Expression scalar2_gluon_vertex =
            new Expression("T^{mn}[k^n,p^n] = -(I/2)*G*(k^n*p^m+k^m*p^n-g^{mn}*(k_a*p^a-m*m))");

    public FeynmanRules() {
        propagator.eval(P.asSubstitution());
    }
}
