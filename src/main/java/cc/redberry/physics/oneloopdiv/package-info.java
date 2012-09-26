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
/**
 * This package contains a number of classes which implement the algorithm for
 * calculating the one-loop counterterms of an arbitrary Lagrangian. The
 * technique of the one-loop background calculations in the general field theory
 * was developed by Konstantin Stepanyantz and can be found in published papers
 * [1-3]. A number of examples are placed in {@code OneLoopCountertermsTest}.
 * Performance benchmarks are placed in {@link  cc.redberry.physics.oneloopdiv.Benchmarks}.
 * The detailed description of the package usage are published [4].
 *
 *
 * <h3>Literature</h3>
 *
 * [1] Petr I. Pronin, Konstantin V. Stepanyantz, <i>One loop counterterms for
 * the dimensional regularization of arbitrary Lagrangians</i>, Phys. Lett. B414
 * (1997) 117-122, <a
 * href="http://arxiv.org/abs/hep-th/9707008">hep-th/9707008</a><br/>
 *
 * [2] P.I. Pronin, K. Stepanyantz, <i>One loop counterterms for higher
 * derivative regularized Lagrangians</i>, Nucl. Phys. B485 (1997) 517-544, <a
 * href="http://arxiv.org/abs/hep-th/9605206">hep-th/9605206</a><br/>
 *
 *
 * [3] Petr I. Pronin, Konstantin V. Stepanyantz, <i>One loop background
 * calculations in the general field theory</i>, <a
 * href="http://arxiv.org/abs/hep-th/9604038">hep-th/9604038</a><br/>
 *
 * [4] D.A. Bolotin and S.V. Poslavsky, <i>Calculation of the one-loop
 * counterterms in an arbitrary theory with Redberry</i>, preparing for
 * publication
 *
 * @see cc.redberry.physics.oneloopdiv.OneLoopCounterterms
 * @see cc.redberry.physics.oneloopdiv.OneLoopInput
 *
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
package cc.redberry.physics.oneloopdiv;
