/*
 * Redberry: symbolic tensor computations.
 *
 * Copyright (c) 2010-2013:
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

package cc.redberry.physics.oneloopdiv;

import cc.redberry.core.indices.IndexType;
import cc.redberry.core.indices.IndicesFactory;
import cc.redberry.core.indices.StructureOfIndices;
import cc.redberry.core.indices.IndicesUtils;
import cc.redberry.core.parser.ParseTokenSimpleTensor;
import cc.redberry.core.parser.preprocessor.IndicesInsertion;
import cc.redberry.core.tensor.Expression;
import cc.redberry.core.tensor.Tensor;
import cc.redberry.core.tensor.Tensors;
import cc.redberry.core.tensor.iterator.TraverseState;
import cc.redberry.core.transformations.EliminateMetricsTransformation;
import cc.redberry.core.transformations.expand.ExpandTransformation;
import cc.redberry.core.transformations.Transformation;
import cc.redberry.core.transformations.Transformer;
import cc.redberry.core.utils.ArraysUtils;
import cc.redberry.core.utils.Indicator;

/**
 * This class is a container of the calculated one-loop counterterms.
 * It has no constructors and can be created using the static
 * method {@link #calculateOneLoopCounterterms(OneLoopInput)}, which
 * performs the whole calculation of the one-loop counterterms.
 * <p>Here is the example of the calculation of one-loop counterterms
 * of vector field: </p>
 * <pre>
 *      //setting symmetries to tensor P
 *      Tensors.addSymmetry("P_\\mu\\nu", IndexType.GreekLower, false, 1, 0);
 *
 *      //input expressions
 *      Expression KINV = Tensors.parseExpression("KINV_\\alpha^\\beta=d_\\alpha^\\beta+ga*n_\\alpha*n^\\beta");
 *      Expression K = Tensors.parseExpression("K^{\\mu\\nu}_\\alpha^{\\beta}=g^{\\mu\\nu}*d_{\\alpha}^{\\beta}-ga/(2*(1+ga))*(g^{\\mu\\beta}*d_\\alpha^\\nu+g^{\\nu\\beta}*d_\\alpha^\\mu)");
 *      Expression S = Tensors.parseExpression("S^\\rho^\\mu_\\nu=0");
 *      Expression W = Tensors.parseExpression("W^{\\alpha}_{\\beta}=P^{\\alpha}_{\\beta}+ga/(2*(1+ga))*R^\\alpha_\\beta");
 *      //F is equal to Riemann for vector field
 *      Expression F = Tensors.parseExpression("F_\\mu\\nu\\alpha\\beta=R_\\mu\\nu\\alpha\\beta");
 *
 *      //tensors M and N are null, since operator order is 2
 *      OneLoopInput input = new OneLoopInput(2, KINV, K, S, W, null, null, F);
 *
 *      //performing the main calculation
 *      OneLoopCounterterms action = OneLoopCounterterms.calculateOneLoopCounterterms(input);
 *      Tensor counterterms = action.counterterms();
 *      //here some transformations can be performed to simplify counterterms
 *      ...
 *      System.out.println(counterterms);
 * </pre>
 * The above code will produce the counterterms, which after some
 * simplifications can be written in form
 * <pre>
 *     (1/24*ga**2+1/4*ga+1/2)*P_\mu\nu*P^\mu\nu + 1/48*ga**2*P**2 + (1/12*ga**2+1/3*ga)*R_\mu\nu*P^\mu\nu +
 *     +(1/24*ga**2+1/12*ga+1/6)*R*P + (1/24*ga**2+1/12*ga-4/15)*R_\mu\nu*R^\mu\nu + (1/48*ga**2+1/12*ga+7/60)*R**2
 * </pre>
 * The divergent part of the one-loop effective action can be obtained by
 * multiplying the resulting counterterms on factor 1/(16*\pi**2*(d-4)) and
 * integrating over the space volume.
 *
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class OneLoopCounterterms {

    private static final String Flat_ =
            "Flat="
                    + "(1/4)*HATS*HATS*HATS*HATS-HATW*HATS*HATS+(1/2)*HATW*HATW+HATS*HATN-HATM+(L-2)*NABLAS_\\mu*HATW^\\mu"
                    + "-L*NABLAS_\\mu*HATW*HATK^\\mu+(1/3)*((L-1)*NABLAS_\\mu^\\mu*HATS*HATS-L*NABLAS_\\mu*HATK^\\mu*HATS*HATS"
                    + "-(L-1)*NABLAS_\\mu*HATS*HATS^\\mu+L*NABLAS_\\mu*HATS*HATS*HATK^\\mu)-(1/2)*NABLAS_\\mu*NABLAS_\\nu*DELTA^{\\mu\\nu}"
                    + "-(1/4)*(L-1)*(L-2)*NABLAS_\\mu*NABLAS_\\nu^{\\mu\\nu}+(1/2)*L*(L-1)*(1/2)*(NABLAS_\\mu*NABLAS_{\\nu }^{\\nu}"
                    + "+NABLAS_{\\nu }*NABLAS_{\\mu }^{\\nu})*HATK^\\mu";
    private static final String WR_ =
            "WR="
                    + "-(1/2)*Power[L,2]*HATW*HATF_{\\mu\\nu}*Kn^\\mu*HATK^\\nu"
                    + "+(1/3)*L*HATW*HATK^\\alpha*DELTA^{\\mu\\nu}*n_\\sigma*R^\\sigma_{\\mu\\alpha\\nu}"
                    + "+(1/3)*Power[L,2]*(L-1)*HATW*HATK^{\\mu\\nu}*HATK^\\alpha*n_\\sigma*R^\\sigma_{\\mu\\alpha\\nu}"
                    + "-(1/6)*(L-2)*(L-3)*HATW^{\\mu\\nu}*R_{\\mu\\nu}";
    private static final String SR_ = "SR=-(1/6)*Power[L,2]*(L-1)*HATS*NABLAF_{\\mu\\alpha\\nu}*Kn^{\\mu\\nu}*HATK^\\alpha"
            + "+(2/3)*L*HATS*NABLAF_{\\mu\\nu\\alpha}*Kn^\\alpha*DELTA^{\\mu\\nu}"
            + "-(1/12)*(L-1)*(L-2)*(L-3)*HATS^{\\alpha\\mu\\nu}*NABLAR_{\\alpha\\mu\\nu}"
            + "-(1/12)*Power[L,2]*(L-1)*(L-2)*HATS*HATK^{\\mu\\nu\\alpha}*HATK^\\beta*n_\\sigma*NABLAR_\\alpha^\\sigma_{\\mu\\beta\\nu}"
            + "+L*(L-1)*HATS*HATK^{\\mu\\nu}*DELTA^{\\alpha\\beta}*n_\\sigma*((5/12)*NABLAR_\\alpha^\\sigma_{\\nu\\beta\\mu}"
            + "-(1/12)*NABLAR_{\\mu}^\\sigma_{\\alpha\\nu\\beta})"
            + "-(1/2)*L*HATS*HATK^\\beta*DELTA^{\\mu\\nu\\alpha}*n_\\sigma*NABLAR_{\\alpha}^{\\sigma}_{\\mu\\beta\\nu}";
    private static final String SSR_ = "SSR=-(1/2)*L*(L-1)*HATS*HATS^\\mu*HATF_{\\mu\\nu}*HATK^{\\nu}+(1/2)*Power[L,2]*HATS*HATS*HATF_{\\mu\\nu}*Kn^{\\mu}*HATK^\\nu"
            + "+(1/12)*(L-1)*(L-2)*HATS*HATS^{\\mu\\nu}*R_{\\mu\\nu}+(1/3)*L*(L-1)*HATS*HATS^\\mu*HATK^\\nu*R_{\\mu\\nu}"
            + "+(1/6)*HATS*HATS*DELTA^{\\mu\\nu}*R_{\\mu\\nu}-(1/6)*L*(L-1)*(L-2)*HATS*HATS^{\\mu\\nu}*HATK^\\alpha*n_\\sigma*R^\\sigma_{\\mu\\alpha\\nu}"
            + "+(1/3)*(L-1)*HATS*HATS^\\alpha*DELTA^{\\mu\\nu}*n_\\sigma*R^\\sigma_{\\mu\\alpha\\nu}"
            + "-(1/3)*Power[L,2]*(L-1)*HATS*HATS*HATK^{\\mu\\nu}*HATK^\\alpha*n_\\sigma*R^\\sigma_{\\mu\\alpha\\nu}"
            + "-(1/3)*L*HATS*HATS*HATK^\\alpha*DELTA^{\\mu\\nu}*n_\\sigma*R^\\sigma_{\\mu\\alpha\\nu}";
    //    public static final String FF_ =
//            "FF="
//            + "-(1/24)*L*L*(L-1)*(L-1)*HATK^{\\mu\\nu}*F_{\\mu\\alpha}*HATK^{\\alpha\\beta}*F_{\\nu\\beta}"
//            + "+(1/24)*L*L*HATK^\\mu*F_{\\beta\\nu}*DELTA^{\\alpha\\beta}*HATK^\\nu*F_{\\alpha\\mu}"
//            + "+(5/24)*L*L*HATK^\\mu*F_{\\beta\\mu}*DELTA^{\\alpha\\beta}*HATK^\\nu*F_{\\alpha\\nu}"
//            + "-(1/48)*L*L*(L-1)*HATK^\\mu*F_{\\beta\\nu}*DELTA^\\nu*HATK^{\\alpha\\beta}*F_{\\alpha\\mu}"
//            + "-(1/48)*L*L*(L-1)*HATK^\\mu*F_{\\beta\\mu}*DELTA^\\nu*HATK^{\\alpha\\beta}*F_{\\alpha\\nu}";
    private static final String FF_ =
            "FF="
                    + "-(1/24)*L*L*(L-1)*(L-1)*HATK^{\\mu\\nu}*F_{\\mu\\alpha}*HATK^{\\alpha\\beta}*F_{\\nu\\beta}"
                    + "+(1/24)*L*L*HATK^\\mu*F_{\\beta\\nu}*DELTA^{\\alpha\\beta}*HATK^\\nu*F_{\\alpha\\mu}"
                    + "-(5/24)*L*L*HATK^\\mu*F_{\\beta\\mu}*DELTA^{\\alpha\\beta}*HATK^\\nu*F_{\\alpha\\nu}"
                    + "-(1/48)*L*L*(L-1)*HATK^\\mu*F_{\\beta\\nu}*DELTA^\\nu*HATK^{\\alpha\\beta}*F_{\\alpha\\mu}"
                    + "-(1/48)*L*L*(L-1)*HATK^\\mu*F_{\\beta\\mu}*DELTA^\\nu*HATK^{\\alpha\\beta}*F_{\\alpha\\nu}";
    private static final String FR_ =
            "FR="
                    + "(1/40)*Power[L,2]*(L-1)*(L-2)*DELTA^\\mu*HATK^\\nu*HATK^{\\alpha\\beta\\gamma}*F_{\\mu\\alpha}*n_\\sigma*R^\\sigma_{\\gamma\\beta\\nu}"
                    + "-Power[L,2]*(L-1)*(L-2)*DELTA^\\nu*HATK^{\\alpha\\beta\\gamma}*HATK^\\mu*n_\\sigma*((1/60)*R^\\sigma_{\\beta\\gamma\\mu}*F_{\\alpha\\nu}"
                    + "+(1/12)*R^\\sigma_{\\beta\\gamma\\nu}*F_{\\alpha\\mu})"
                    + "+Power[L,2]*Power[(L-1),2]*DELTA^\\alpha*HATK^{\\beta\\gamma}*HATK^{\\mu\\nu}*n_\\sigma*((1/60)*R^\\sigma_{\\beta\\mu\\gamma}*F_{\\alpha\\nu}"
                    + "+(1/20)*R^\\sigma_{\\alpha\\mu\\gamma}*F_{\\nu\\beta}+(1/15)*R^\\sigma_{\\gamma\\mu\\alpha}*F_{\\nu\\beta}"
                    + "+(1/60)*R^\\sigma_{\\mu\\nu\\gamma}*F_{\\alpha\\beta})+Power[L,2]*(L-1)*DELTA^{\\alpha\\beta}*HATK^{\\gamma\\delta}*HATK^{\\mu}"
                    + "*n_\\sigma*((4/15)*R^\\sigma_{\\delta\\beta\\gamma}*F_{\\alpha\\mu}-(1/30)*R^\\sigma_{\\beta\\delta\\alpha}*F_{\\gamma\\mu}"
                    + "-(1/15)*R^\\sigma_{\\alpha\\gamma\\mu}*F_{\\beta\\delta}-(1/30)*R^\\sigma_{\\gamma\\alpha\\mu}*F_{\\beta\\delta})"
                    + "+Power[L,2]*(L-1)*DELTA^{\\alpha\\beta}*HATK^\\gamma*HATK^{\\mu\\nu}*n_\\sigma*((7/60)*R^\\sigma_{\\alpha\\beta\\mu}*F_{\\gamma\\nu}"
                    + "-(11/60)*R^\\sigma_{\\beta\\mu\\gamma}*F_{\\alpha\\nu}+(1/5)*R^\\sigma_{\\mu\\alpha\\gamma}*F_{\\beta\\nu}"
                    + "+(1/60)*R^\\sigma_{\\mu\\alpha\\nu}*F_{\\gamma\\beta})"
                    + "+Power[L,2]*DELTA^{\\mu\\alpha\\beta}*HATK^\\gamma*HATK^\\nu*n_\\sigma"
                    + "*((7/20)*R^\\sigma_{\\alpha\\gamma\\beta}*F_{\\nu\\mu}+(1/10)*R^\\sigma_{\\alpha\\beta\\nu}*F_{\\gamma\\mu})";
    //    public static final String RR_ =
//            "RR="
//          + "+Power[L,2]*(L-1)*HATK^\\delta*DELTA^{\\alpha\\beta\\gamma}*HATK^{\\mu\\nu}*n_\\sigma*n_\\rho*"
//            + "((-1/10)*R^\\rho_{\\mu\\gamma\\nu}*R^\\sigma_{\\alpha\\delta\\beta}+(1/15)*R^\\rho_{\\delta\\alpha\\nu}*R^\\sigma_{\\beta\\mu\\gamma}+(1/60)*R^\\rho_{\\beta\\delta\\nu}*R^\\sigma_{\\gamma\\mu\\alpha})"
//           
//            ;
    private static final String RR_ =
            "RR="
                    + "(1/10)*Power[L,2]*HATK^\\delta*DELTA^{\\mu\\nu\\alpha\\beta}*HATK^\\gamma*n_\\sigma*n_\\rho*R^\\sigma_{\\alpha\\beta\\gamma}*R^\\rho_{\\mu\\nu\\delta}"
                    + "+Power[L,2]*Power[(L-1),2]*(L-2)*HATK^{\\beta\\gamma\\delta}*DELTA^\\alpha*HATK^{\\mu\\nu}*n_\\sigma*n_\\rho*"
                    + "((2/45)*R^\\rho_{\\alpha\\delta\\nu}*R^\\sigma_{\\beta\\mu\\gamma}-(1/120)*R^\\rho_{\\delta\\alpha\\nu}*R^\\sigma_{\\beta\\mu\\gamma})"
                    + "+Power[L,2]*(L-1)*HATK^\\delta*DELTA^{\\alpha\\beta\\gamma}*HATK^{\\mu\\nu}*n_\\sigma*n_\\rho*"
                    + "((-1/10)*R^\\rho_{\\mu\\gamma\\nu}*R^\\sigma_{\\alpha\\delta\\beta}+(1/15)*R^\\rho_{\\delta\\alpha\\nu}*R^\\sigma_{\\beta\\mu\\gamma}+(1/60)*R^\\rho_{\\beta\\delta\\nu}*R^\\sigma_{\\gamma\\mu\\alpha})"
                    + "+Power[L,2]*Power[(L-1),2]*HATK^{\\gamma\\delta}*DELTA^{\\alpha\\beta}*HATK^{\\mu\\nu}*n_\\sigma*n_\\rho*"
                    + "(-(1/20)*R^\\rho_{\\mu\\beta\\nu}*R^\\sigma_{\\delta\\alpha\\gamma}+(1/180)*R^\\rho_{\\alpha\\nu\\beta}*R^\\sigma_{\\gamma\\delta\\mu}-(7/360)*R^\\rho_{\\mu\\gamma\\nu}*R^\\sigma_{\\alpha\\delta\\beta}-(1/240)*R^\\rho_{\\delta\\beta\\nu}*R^\\sigma_{\\gamma\\alpha\\mu}-(1/120)*R^\\rho_{\\beta\\gamma\\nu}*R^\\sigma_{\\alpha\\delta\\mu}-(1/30)*R^\\rho_{\\delta\\beta\\nu}*R^\\sigma_{\\alpha\\gamma\\mu})"
                    + "+Power[L,2]*(L-1)*(L-2)*HATK^\\delta*DELTA^{\\mu\\nu}*HATK^{\\alpha\\beta\\gamma}*n_\\sigma*n_\\rho*"
                    + "((-1/30)*R^\\rho_{\\gamma\\nu\\beta}*R^\\sigma_{\\alpha\\delta\\mu}-(1/180)*R^\\rho_{\\mu\\gamma\\nu}*R^\\sigma_{\\alpha\\beta\\delta}+(1/180)*R^\\rho_{\\mu\\gamma\\delta}*R^\\sigma_{\\alpha\\beta\\nu})"
                    + "+Power[L,2]*Power[(L-1),2]*(L-2)*HATK^{\\mu\\nu}*DELTA^{\\delta}*HATK^{\\alpha\\beta\\gamma}*n_\\sigma*n_\\rho*"
                    + "((1/45)*R^\\rho_{\\mu\\gamma\\nu}*R^\\sigma_{\\alpha\\beta\\delta}-(1/80)*R^\\rho_{\\beta\\nu\\gamma}*R^\\sigma_{\\mu\\alpha\\delta}+(1/90)*R^\\rho_{\\beta\\nu\\gamma}*R^\\sigma_{\\delta\\alpha\\mu})"
                    + "+Power[L,2]*(L-1)*HATK^{\\mu\\nu}*DELTA^{\\alpha\\beta\\gamma}*HATK^\\delta*n_\\sigma*n_\\rho*"
                    + "((7/120)*R^\\rho_{\\beta\\gamma\\nu}*R^\\sigma_{\\mu\\alpha\\delta}-(3/40)*R^\\rho_{\\beta\\gamma\\delta}*R^\\sigma_{\\mu\\alpha\\nu}+(1/120)*R^\\rho_{\\delta\\gamma\\nu}*R^\\sigma_{\\alpha\\beta\\mu})"
                    + "+Power[L,2]*(L-1)*(L-2)*HATK^{\\alpha\\beta\\gamma}*DELTA^{\\mu\\nu}*HATK^\\delta*n_\\sigma*n_\\rho*"
                    + "(-(1/24)*R^\\rho_{\\mu\\gamma\\nu}*R^\\sigma_{\\alpha\\beta\\delta}-(1/180)*R^\\rho_{\\nu\\gamma\\delta}*R^\\sigma_{\\alpha\\beta\\mu}-(1/360)*R^\\rho_{\\delta\\gamma\\nu}*R^\\sigma_{\\alpha\\beta\\mu})"
                    + "-(1/120)*Power[L,2]*(L-1)*(L-2)*(L-3)*HATK^{\\mu\\nu\\alpha\\beta}*DELTA^{\\delta}*HATK^\\gamma*n_\\sigma*n_\\rho*R^\\rho_{\\alpha\\beta\\gamma}*R^\\sigma_{\\mu\\nu\\delta}"
                    + "-(1/80)*Power[L,2]*Power[(L-1),2]*(L-2)*(L-3)*HATK^{\\alpha\\beta\\gamma\\delta}*HATK^{\\mu\\nu}*n_\\sigma*n_\\rho*R^\\rho_{\\beta\\gamma\\mu}*R^\\sigma_{\\alpha\\delta\\nu}"
                    + "+Power[L,2]*HATK^\\mu*DELTA^{\\alpha\\beta\\gamma}*HATK^\\nu*n_\\rho*(-(1/8)*R_{\\beta\\gamma}*R^\\rho_{\\nu\\alpha\\mu}+(3/20)*R_{\\beta\\gamma}*R^\\rho_{\\mu\\alpha\\nu}+(3/40)*R_{\\alpha\\mu}*R^\\rho_{\\beta\\gamma\\nu}+(1/40)*R^\\sigma_{\\beta\\gamma\\mu}*R^\\rho_{\\nu\\alpha\\sigma}-(3/20)*R^\\sigma_{\\alpha\\beta\\mu}*R^\\rho_{\\gamma\\nu\\sigma}+(1/10)*R^\\sigma_{\\alpha\\beta\\nu}*R^\\rho_{\\gamma\\mu\\sigma})"
                    + "+Power[L,2]*(L-1)*HATK^\\gamma*DELTA^{\\alpha\\beta}*HATK^{\\mu\\nu}*n_\\rho*"
                    + "((1/20)*R_{\\alpha\\nu}*R^\\rho_{\\gamma\\beta\\mu}+(1/20)*R_{\\alpha\\gamma}*R^\\rho_{\\mu\\beta\\nu}+(1/10)*R_{\\alpha\\beta}*R^\\rho_{\\mu\\gamma\\nu}+(1/20)*R^\\sigma_{\\alpha\\nu\\gamma}*R^\\rho_{\\sigma\\beta\\mu}-(1/60)*R^\\sigma_{\\mu\\alpha\\nu}*R^\\rho_{\\beta\\sigma\\gamma}+(1/10)*R^\\sigma_{\\alpha\\beta\\gamma}*R^\\rho_{\\mu\\sigma\\nu}-(1/12)*R^\\sigma_{\\alpha\\beta\\nu}*R^\\rho_{\\mu\\sigma\\gamma})"
                    + "+Power[L,2]*Power[(L-1),2]*HATK^{\\alpha\\beta}*DELTA^{\\gamma}*HATK^{\\mu\\nu}*n_\\rho*"
                    + "((1/60)*R_{\\alpha\\mu}*R^\\rho_{\\beta\\nu\\gamma}-(1/20)*R_{\\alpha\\mu}*R^\\rho_{\\gamma\\nu\\beta}+(1/120)*R_{\\alpha\\beta}*R^\\rho_{\\mu\\nu\\gamma}+(3/40)*R_{\\alpha\\gamma}*R^\\rho_{\\nu\\beta\\mu}+(1/20)*R^\\sigma_{\\gamma\\mu\\alpha}*R^\\rho_{\\nu\\sigma\\beta}+(1/120)*R^\\sigma_{\\alpha\\mu\\gamma}*R^\\rho_{\\beta\\nu\\sigma}-(1/40)*R^\\sigma_{\\alpha\\mu\\gamma}*R^\\rho_{\\sigma\\nu\\beta}+(1/40)*R^\\sigma_{\\alpha\\mu\\beta}*R^\\rho_{\\sigma\\nu\\gamma}-(1/20)*R^\\sigma_{\\alpha\\mu\\beta}*R^\\rho_{\\gamma\\nu\\sigma}-(1/40)*R^\\sigma_{\\mu\\beta\\nu}*R^\\rho_{\\gamma\\sigma\\alpha})"
                    + "+Power[L,2]*(L-1)*HATK^{\\alpha\\beta}*DELTA^{\\mu\\nu}*HATK^{\\gamma}*n_\\rho*"
                    + "((1/20)*R^\\sigma_{\\mu\\nu\\beta}*R^\\rho_{\\gamma\\sigma\\alpha}-(7/60)*R^\\sigma_{\\beta\\mu\\alpha}*R^\\rho_{\\gamma\\nu\\sigma}+(1/20)*R^\\sigma_{\\beta\\mu\\alpha}*R^\\rho_{\\sigma\\nu\\gamma}+(1/10)*R^\\sigma_{\\mu\\beta\\gamma}*R^\\rho_{\\nu\\alpha\\sigma}+(1/60)*R^\\sigma_{\\beta\\mu\\gamma}*R^\\rho_{\\alpha\\nu\\sigma}+(7/120)*R_{\\alpha\\beta}*R^\\rho_{\\nu\\gamma\\mu}+(11/60)*R_{\\beta\\mu}*R^\\rho_{\\nu\\alpha\\gamma})"
                    + "+Power[L,2]*(L-1)*(L-2)*HATK^{\\alpha\\beta\\gamma}*DELTA^{\\mu}*HATK^{\\nu}*n_\\rho*"
                    + "((7/240)*R_{\\alpha\\beta}*R^\\rho_{\\gamma\\mu\\nu}+(7/240)*R_{\\alpha\\nu}*R^\\rho_{\\beta\\gamma\\mu}-(1/60)*R_{\\alpha\\mu}*R^\\rho_{\\beta\\gamma\\nu}-(1/24)*R^\\sigma_{\\alpha\\beta\\nu}*R^\\rho_{\\sigma\\gamma\\mu}+(1/15)*R^\\sigma_{\\alpha\\beta\\nu}*R^\\rho_{\\mu\\gamma\\sigma}+(1/40)*R^\\sigma_{\\alpha\\beta\\mu}*R^\\rho_{\\sigma\\gamma\\nu}+(1/40)*R_{\\beta\\gamma}*R^\\rho_{\\nu\\mu\\alpha}+(1/48)*R^\\sigma_{\\beta\\gamma\\mu}*R^\\rho_{\\nu\\alpha\\sigma})"
                    + "+Power[L,2]*Power[(L-1),2]*(L-2)*HATK^{\\alpha\\beta\\gamma}*HATK^{\\mu\\nu}*n_\\rho*"
                    + "((-7/240)*R_{\\alpha\\mu}*R^\\rho_{\\beta\\gamma\\nu}+(1/240)*R_{\\beta\\gamma}*R^\\rho_{\\mu\\alpha\\nu}-(1/40)*R^\\sigma_{\\alpha\\mu\\beta}*R^\\rho_{\\nu\\gamma\\sigma})"
                    + "+L*(L-1)*(L-2)*(L-3)*HATK^{\\mu\\nu\\alpha\\beta}*"
                    + "((1/180)*R_{\\mu\\nu}*R_{\\alpha\\beta}+(7/720)*R^\\sigma_{\\alpha\\beta\\rho}*R^\\rho_{\\mu\\nu\\sigma})";
    //    public static final String RR_ =
//            "RR="
//            + "(1/10)*Power[L,2]*HATK^\\delta*DELTA^{\\mu\\nu\\alpha\\beta}*HATK^\\gamma*n_\\sigma*n_\\rho*R^\\sigma_{\\alpha\\beta\\gamma}*R^\\rho_{\\mu\\nu\\delta}"
//            /*
//             * (L-2)
//             */ + "+Power[L,2]*Power[(L-1),2]*(L-2)*HATK^{\\beta\\gamma\\delta}*DELTA^\\alpha*HATK^{\\mu\\nu}*n_\\sigma*n_\\rho*"
//            + "((2/45)*R^\\rho_{\\alpha\\delta\\nu}*R^\\sigma_{\\beta\\mu\\gamma}-(1/120)*R^\\rho_{\\delta\\alpha\\nu}*R^\\sigma_{\\beta\\mu\\gamma})"
//            + "+Power[L,2]*(L-1)*HATK^\\delta*DELTA^{\\alpha\\beta\\gamma}*HATK^{\\mu\\nu}*n_\\sigma*n_\\rho*"
//            + "((-1/10)*R^\\rho_{\\mu\\gamma\\nu}*R^\\sigma_{\\alpha\\delta\\beta}+(1/15)*R^\\rho_{\\delta\\alpha\\nu}*R^\\sigma_{\\beta\\mu\\gamma}+(1/60)*R^\\rho_{\\beta\\delta\\nu}*R^\\sigma_{\\gamma\\mu\\alpha})"
//            + "+Power[L,2]*Power[(L-1),2]*HATK^{\\gamma\\delta}*DELTA^{\\alpha\\beta}*HATK^{\\mu\\nu}*n_\\sigma*n_\\rho*"
//            + "(-(1/20)*R^\\rho_{\\mu\\beta\\nu}*R^\\sigma_{\\delta\\alpha\\gamma}+(1/180)*R^\\rho_{\\alpha\\nu\\beta}*R^\\sigma_{\\gamma\\delta\\mu}-(7/360)*R^\\rho_{\\mu\\gamma\\nu}*R^\\sigma_{\\alpha\\delta\\beta}-(1/240)*R^\\rho_{\\delta\\beta\\nu}*R^\\sigma_{\\gamma\\alpha\\mu}-(1/120)*R^\\rho_{\\beta\\gamma\\nu}*R^\\sigma_{\\alpha\\delta\\mu}-(1/30)*R^\\rho_{\\delta\\beta\\nu}*R^\\sigma_{\\alpha\\gamma\\mu})"
//            /*
//             * (L-2)
//             */ + "+Power[L,2]*(L-1)*(L-2)*HATK^\\delta*DELTA^{\\mu\\nu}*HATK^{\\alpha\\beta\\gamma}*n_\\sigma*n_\\rho*"
//            + "((-1/30)*R^\\rho_{\\gamma\\nu\\beta}*R^\\sigma_{\\alpha\\delta\\mu}-(1/180)*R^\\rho_{\\mu\\gamma\\nu}*R^\\sigma_{\\alpha\\beta\\delta}+(1/180)*R^\\rho_{\\mu\\gamma\\delta}*R^\\sigma_{\\alpha\\beta\\nu})"
//            /*
//             * (L-2)
//             */ + "+Power[L,2]*Power[(L-1),2]*(L-2)*HATK^{\\mu\\nu}*DELTA^{\\delta}*HATK^{\\alpha\\beta\\gamma}*n_\\sigma*n_\\rho*"
//            + "((1/45)*R^\\rho_{\\mu\\gamma\\nu}*R^\\sigma_{\\alpha\\beta\\delta}-(1/80)*R^\\rho_{\\beta\\nu\\gamma}*R^\\sigma_{\\mu\\alpha\\delta}+(1/90)*R^\\rho_{\\beta\\nu\\gamma}*R^\\sigma_{\\delta\\alpha\\mu})"
//            + "+Power[L,2]*(L-1)*HATK^{\\mu\\nu}*DELTA^{\\alpha\\beta\\gamma}*HATK^\\delta*n_\\sigma*n_\\rho*"
//            + "((7/120)*R^\\rho_{\\beta\\gamma\\nu}*R^\\sigma_{\\mu\\alpha\\delta}-(3/40)*R^\\rho_{\\beta\\gamma\\delta}*R^\\sigma_{\\mu\\alpha\\nu}+(1/120)*R^\\rho_{\\delta\\gamma\\nu}*R^\\sigma_{\\alpha\\beta\\mu})"
//            /*
//             * (L-2)
//             */ + "+Power[L,2]*(L-1)*(L-2)*HATK^{\\alpha\\beta\\gamma}*DELTA^{\\mu\\nu}*HATK^\\delta*n_\\sigma*n_\\rho*"
//            + "(-(1/24)*R^\\rho_{\\mu\\gamma\\nu}*R^\\sigma_{\\alpha\\beta\\delta}-(1/180)*R^\\rho_{\\nu\\gamma\\delta}*R^\\sigma_{\\alpha\\beta\\mu}-(1/360)*R^\\rho_{\\delta\\gamma\\nu}*R^\\sigma_{\\alpha\\beta\\mu})"
//            /*
//             * (L-2)
//             *//*
//             * kv reduce
//             */ + "-(1/120)*Power[L,2]*(L-1)*(L-2)*(L-3)*HATK^\\delta*DELTA^{\\gamma}*HATK^{\\mu\\nu\\alpha\\beta}*n_\\sigma*n_\\rho*R^\\rho_{\\alpha\\beta\\gamma}*R^\\sigma_{\\mu\\nu\\delta}"
//            ///*(L-2)*//*kv paper*/    + "-(1/120)*Power[L,2]*(L-1)*(L-2)*(L-3)*HATK^{\\mu\\nu\\alpha\\beta}*DELTA^{\\delta}*HATK^\\gamma*n_\\sigma*n_\\rho*R^\\rho_{\\alpha\\beta\\gamma}*R^\\sigma_{\\mu\\nu\\delta}"
//
//            /*
//             * (L-2)
//             *//*
//             * kv reduce
//             */ + "-(1/80)*Power[L,2]*Power[(L-1),2]*(L-2)*(L-3)*HATK^{\\alpha\\beta\\gamma\\delta}*HATK^{\\mu\\nu}*n_\\sigma*n_\\rho*R^\\rho_{\\gamma\\beta\\mu}*R^\\sigma_{\\alpha\\delta\\nu}"
//            ///*(L-2)*//*kv paper*/    + "-(1/80)*Power[L,2]*Power[(L-1),2]*(L-2)*(L-3)*HATK^{\\alpha\\beta\\gamma\\delta}*HATK^{\\mu\\nu}*n_\\sigma*n_\\rho*R^\\rho_{\\beta\\gamma\\mu}*R^\\sigma_{\\alpha\\delta\\nu}"
//
//            + "+Power[L,2]*HATK^\\mu*DELTA^{\\alpha\\beta\\gamma}*HATK^\\nu*n_\\rho*(-(1/8)*R_{\\beta\\gamma}*R^\\rho_{\\nu\\alpha\\mu}+(3/20)*R_{\\beta\\gamma}*R^\\rho_{\\mu\\alpha\\nu}+(3/40)*R_{\\alpha\\mu}*R^\\rho_{\\beta\\gamma\\nu}+(1/40)*R^\\sigma_{\\beta\\gamma\\mu}*R^\\rho_{\\nu\\alpha\\sigma}-(3/20)*R^\\sigma_{\\alpha\\beta\\mu}*R^\\rho_{\\gamma\\nu\\sigma}+(1/10)*R^\\sigma_{\\alpha\\beta\\nu}*R^\\rho_{\\gamma\\mu\\sigma})"
//            + "+Power[L,2]*(L-1)*HATK^\\gamma*DELTA^{\\alpha\\beta}*HATK^{\\mu\\nu}*n_\\rho*"
//            + "((1/20)*R_{\\alpha\\nu}*R^\\rho_{\\gamma\\beta\\mu}+(1/20)*R_{\\alpha\\gamma}*R^\\rho_{\\mu\\beta\\nu}+(1/10)*R_{\\alpha\\beta}*R^\\rho_{\\mu\\gamma\\nu}+(1/20)*R^\\sigma_{\\alpha\\nu\\gamma}*R^\\rho_{\\sigma\\beta\\mu}-(1/60)*R^\\sigma_{\\mu\\alpha\\nu}*R^\\rho_{\\beta\\sigma\\gamma}+(1/10)*R^\\sigma_{\\alpha\\beta\\gamma}*R^\\rho_{\\mu\\sigma\\nu}-(1/12)*R^\\sigma_{\\alpha\\beta\\nu}*R^\\rho_{\\mu\\sigma\\gamma})"
//            + "+Power[L,2]*Power[(L-1),2]*HATK^{\\alpha\\beta}*DELTA^{\\gamma}*HATK^{\\mu\\nu}*n_\\rho*"
//            + "((1/60)*R_{\\alpha\\mu}*R^\\rho_{\\beta\\nu\\gamma}-(1/20)*R_{\\alpha\\mu}*R^\\rho_{\\gamma\\nu\\beta}+(1/120)*R_{\\alpha\\beta}*R^\\rho_{\\mu\\nu\\gamma}+(3/40)*R_{\\alpha\\gamma}*R^\\rho_{\\nu\\beta\\mu}+(1/20)*R^\\sigma_{\\gamma\\mu\\alpha}*R^\\rho_{\\nu\\sigma\\beta}+(1/120)*R^\\sigma_{\\alpha\\mu\\gamma}*R^\\rho_{\\beta\\nu\\sigma}-(1/40)*R^\\sigma_{\\alpha\\mu\\gamma}*R^\\rho_{\\sigma\\nu\\beta}+(1/40)*R^\\sigma_{\\alpha\\mu\\beta}*R^\\rho_{\\sigma\\nu\\gamma}-(1/20)*R^\\sigma_{\\alpha\\mu\\beta}*R^\\rho_{\\gamma\\nu\\sigma}-(1/40)*R^\\sigma_{\\mu\\beta\\nu}*R^\\rho_{\\gamma\\sigma\\alpha})"
//            + "+Power[L,2]*(L-1)*HATK^{\\alpha\\beta}*DELTA^{\\mu\\nu}*HATK^{\\gamma}*n_\\rho*"
//            + "((1/20)*R^\\sigma_{\\mu\\nu\\beta}*R^\\rho_{\\gamma\\sigma\\alpha}-(7/60)*R^\\sigma_{\\beta\\mu\\alpha}*R^\\rho_{\\gamma\\nu\\sigma}+(1/20)*R^\\sigma_{\\beta\\mu\\alpha}*R^\\rho_{\\sigma\\nu\\gamma}+(1/10)*R^\\sigma_{\\mu\\beta\\gamma}*R^\\rho_{\\nu\\alpha\\sigma}+(1/60)*R^\\sigma_{\\beta\\mu\\gamma}*R^\\rho_{\\alpha\\nu\\sigma}+(7/120)*R_{\\alpha\\beta}*R^\\rho_{\\nu\\gamma\\mu}+(11/60)*R_{\\beta\\mu}*R^\\rho_{\\nu\\alpha\\gamma})"
//            /*
//             * (L-2)
//             */ + "+Power[L,2]*(L-1)*(L-2)*HATK^{\\alpha\\beta\\gamma}*DELTA^{\\mu}*HATK^{\\nu}*n_\\rho*"
//            + "((7/240)*R_{\\alpha\\beta}*R^\\rho_{\\gamma\\mu\\nu}+(7/240)*R_{\\alpha\\nu}*R^\\rho_{\\beta\\gamma\\mu}-(1/60)*R_{\\alpha\\mu}*R^\\rho_{\\beta\\gamma\\nu}-(1/24)*R^\\sigma_{\\alpha\\beta\\nu}*R^\\rho_{\\sigma\\gamma\\mu}+(1/15)*R^\\sigma_{\\alpha\\beta\\nu}*R^\\rho_{\\mu\\gamma\\sigma}+(1/40)*R^\\sigma_{\\alpha\\beta\\mu}*R^\\rho_{\\sigma\\gamma\\nu}+(1/40)*R_{\\beta\\gamma}*R^\\rho_{\\nu\\mu\\alpha}+(1/48)*R^\\sigma_{\\beta\\gamma\\mu}*R^\\rho_{\\nu\\alpha\\sigma})"
//            /*
//             * (L-2)
//             *//*
//             * kv reduce
//             */ + "+Power[L,2]*Power[(L-1),2]*(L-2)*HATK^{\\mu\\nu}*HATK^{\\alpha\\beta\\gamma}*n_\\rho*"
//            ///*(L-2)*//*kv paper*/   + "+Power[L,2]*Power[(L-1),2]*(L-2)*HATK^{\\alpha\\beta\\gamma}*HATK^{\\mu\\nu}*n_\\rho*"
//            + "((-7/240)*R_{\\alpha\\mu}*R^\\rho_{\\beta\\gamma\\nu}+(1/240)*R_{\\beta\\gamma}*R^\\rho_{\\mu\\alpha\\nu}-(1/40)*R^\\sigma_{\\alpha\\mu\\beta}*R^\\rho_{\\nu\\gamma\\sigma})"
//            /*
//             * (L-2)
//             */ + "+L*(L-1)*(L-2)*(L-3)*HATK^{\\mu\\nu\\alpha\\beta}*"
//            + "((1/180)*R_{\\mu\\nu}*R_{\\alpha\\beta}+(7/720)*R^\\sigma_{\\alpha\\beta\\rho}*R^\\rho_{\\mu\\nu\\sigma})";
    //From KV
    //public static final String RR_ = "RR=L**2/10 *(R^{\\sigma}_{\\alpha \\beta \\gamma }*R^{\\rho}_{\\mu \\nu \\delta } *n_{\\sigma}*n_{\\rho})*HATK^{\\delta  }*DELTA^{\\mu \\nu \\alpha \\beta  }*HATK^{\\gamma  } + L**2*(L-1)**2*(L-2)*n_{\\sigma}*n_{\\rho} *(2/45*R^{\\rho}_{\\alpha \\delta \\nu }*R^{\\sigma}_{\\beta \\mu \\gamma }-1/120*R^{\\rho}_{\\delta \\alpha \\nu }*R^{\\sigma}_{\\beta \\mu \\gamma }) *HATK^{\\beta \\gamma \\delta  }*DELTA^{\\alpha  }*HATK^{\\mu \\nu  } + L**2*(L-1)*n_{\\rho}*n_{\\sigma} *(-1/10*R^{\\sigma}_{\\mu \\gamma \\nu }*R^{\\rho}_{\\alpha \\delta \\beta }+1/15*R^{\\sigma}_{\\delta \\alpha \\nu }*R^{\\rho}_{\\beta \\mu \\gamma } +1/60*R^{\\sigma}_{\\beta \\delta \\nu }*R^{\\rho}_{\\gamma \\mu \\alpha }) *HATK^{\\delta  }*DELTA^{\\alpha \\beta \\gamma  }*HATK^{\\mu \\nu  } + L**2*(L-1)**2*n_{\\sigma}*n_{\\rho} *(-1/20*R^{\\rho}_{\\mu \\beta \\nu }*R^{\\sigma}_{\\delta \\alpha \\gamma }+1/180*R^{\\rho}_{\\alpha \\nu \\beta }*R^{\\sigma}_{\\gamma \\delta \\mu } -7/360*R^{\\rho}_{\\mu \\gamma \\nu }*R^{\\sigma}_{\\alpha \\delta \\beta }-1/240*R^{\\rho}_{\\delta \\beta \\nu }*R^{\\sigma}_{\\gamma \\alpha \\mu } -1/120*R^{\\rho}_{\\beta \\gamma \\nu }*R^{\\sigma}_{\\alpha \\delta \\mu }-1/30*R^{\\rho}_{\\delta \\beta \\nu }*R^{\\sigma}_{\\alpha \\gamma \\mu }) *HATK^{\\gamma \\delta  }*DELTA^{\\alpha \\beta  }*HATK^{\\mu \\nu  } + L**2*(L-1)*(L-2)*n_{\\sigma}*n_{\\rho} *(-1/30*R^{\\sigma}_{\\gamma \\nu \\beta }*R^{\\rho}_{\\alpha \\delta \\mu }-1/180*R^{\\sigma}_{\\mu \\gamma \\nu }*R^{\\rho}_{\\alpha \\beta \\delta } +1/180*R^{\\sigma}_{\\mu \\gamma \\delta }*R^{\\rho}_{\\alpha \\beta \\nu }) *HATK^{\\delta  }*DELTA^{\\mu \\nu  }*HATK^{\\alpha \\beta \\gamma  } + L**2*(L-1)**2*(L-2)*n_{\\sigma}*n_{\\rho} *(1/45*R^{\\rho}_{\\mu \\gamma \\nu }*R^{\\sigma}_{\\alpha \\beta \\delta }-1/80*R^{\\rho}_{\\beta \\nu \\gamma }*R^{\\sigma}_{\\mu \\alpha \\delta } +1/90*R^{\\rho}_{\\beta \\nu \\gamma }*R^{\\sigma}_{\\delta \\alpha \\mu }) *HATK^{\\mu \\nu  }*DELTA^{\\delta  }*HATK^{\\alpha \\beta \\gamma  } + L**2*(L-1)*n_{\\sigma}*n_{\\rho} *(7/120*R^{\\rho}_{\\beta \\gamma \\nu }*R^{\\sigma}_{\\mu \\alpha \\delta }-3/40*R^{\\rho}_{\\beta \\gamma \\delta }*R^{\\sigma}_{\\mu \\alpha \\nu } +1/120*R^{\\rho}_{\\delta \\gamma \\nu }*R^{\\sigma}_{\\alpha \\beta \\mu }) *HATK^{\\mu \\nu  }*DELTA^{\\alpha \\beta \\gamma  }*HATK^{\\delta  } + L**2*(L-1)*(L-2)*n_{\\sigma}*n_{\\rho} *(-1/24*R^{\\rho}_{\\mu \\gamma \\nu }*R^{\\sigma}_{\\alpha \\beta \\delta }-1/180*R^{\\rho}_{\\nu \\gamma \\delta }*R^{\\sigma}_{\\alpha \\beta \\mu } -1/360*R^{\\rho}_{\\delta \\gamma \\nu }*R^{\\sigma}_{\\alpha \\beta \\mu }) *HATK^{\\alpha \\beta \\gamma  }*DELTA^{\\mu \\nu  }*HATK^{\\delta  } - L**2*(L-1)*(L-2)*(L-3)*(n_{\\sigma}*n_{\\rho} *R^{\\sigma}_{\\alpha \\beta \\gamma }*R^{\\rho}_{\\mu \\nu \\delta }) *HATK^{\\delta  }*DELTA^{\\gamma  }*HATK^{\\mu \\nu \\alpha \\beta  } /120 - L**2*(L-1)**2*(L-2)*(L-3)*(n_{\\sigma}*n_{\\rho} *R^{\\rho}_{\\gamma \\beta \\mu }*R^{\\sigma}_{\\alpha \\delta \\nu }) *HATK^{\\alpha \\beta \\gamma \\delta  }*HATK^{\\mu \\nu  } /80 + L**2*n_{\\rho} *(-1/8*R_{\\beta \\gamma}*R^{\\rho}_{\\nu \\alpha \\mu }+3/20*R_{\\beta \\gamma}*R^{\\rho}_{\\mu \\alpha \\nu } +3/40*R_{\\alpha \\mu}*R^{\\rho}_{\\beta \\gamma \\nu }+1/40*R^{\\sigma}_{\\beta \\gamma \\mu }*R^{\\rho}_{\\nu \\alpha \\sigma } -3/20*R^{\\sigma}_{\\alpha \\beta \\mu }*R^{\\rho}_{\\gamma \\nu \\sigma }+1/10*R^{\\sigma}_{\\alpha \\beta \\nu }*R^{\\rho}_{\\gamma \\mu \\sigma }) *HATK^{\\mu  }*DELTA^{\\alpha \\beta \\gamma  }*HATK^{\\nu  } + L**2*(L-1)*n_{\\rho} *(1/20*R_{\\alpha \\nu}*R^{\\rho}_{\\gamma \\beta \\mu } +1/20*R_{\\alpha \\gamma}*R^{\\rho}_{\\mu \\beta \\nu }+1/10*R_{\\alpha \\beta}*R^{\\rho}_{\\mu \\gamma \\nu } +1/20*R^{\\sigma}_{\\alpha \\nu \\gamma }*R^{\\rho}_{\\sigma \\beta \\mu }-1/60*R^{\\sigma}_{\\mu \\alpha \\nu }*R^{\\rho}_{\\beta \\sigma \\gamma } +1/10*R^{\\sigma}_{\\alpha \\beta \\gamma }*R^{\\rho}_{\\mu \\sigma \\nu }-1/12*R^{\\sigma}_{\\alpha \\beta \\nu }*R^{\\rho}_{\\mu \\sigma \\gamma }) *HATK^{\\gamma  }*DELTA^{\\alpha \\beta  }*HATK^{\\mu \\nu  } + L**2*(L-1)**2*n_{\\rho} *(1/60*R_{\\alpha \\mu}*R^{\\rho}_{\\beta \\nu \\gamma }-1/20*R_{\\alpha \\mu}*R^{\\rho}_{\\gamma \\nu \\beta } +1/120*R_{\\alpha \\beta}*R^{\\rho}_{\\mu \\nu \\gamma }+3/40*R_{\\alpha \\gamma}*R^{\\rho}_{\\nu \\beta \\mu } +1/20*R^{\\sigma}_{\\gamma \\mu \\alpha }*R^{\\rho}_{\\nu \\sigma \\beta }+1/120*R^{\\sigma}_{\\alpha \\mu \\gamma }*R^{\\rho}_{\\beta \\nu \\sigma } -1/40*R^{\\sigma}_{\\alpha \\mu \\gamma }*R^{\\rho}_{\\sigma \\nu \\beta }+1/40*R^{\\sigma}_{\\alpha \\mu \\beta }*R^{\\rho}_{\\sigma \\nu \\gamma } -1/20*R^{\\sigma}_{\\alpha \\mu \\beta }*R^{\\rho}_{\\gamma \\nu \\sigma }-1/40*R^{\\sigma}_{\\mu \\beta \\nu }*R^{\\rho}_{\\gamma \\sigma \\alpha }) *HATK^{\\alpha \\beta  }*DELTA^{\\gamma  }*HATK^{\\mu \\nu  } + L**2*(L-1)*n_{\\rho} *(1/20*R^{\\sigma}_{\\mu \\nu \\beta }*R^{\\rho}_{\\gamma \\sigma \\alpha }-7/60*R^{\\sigma}_{\\beta \\mu \\alpha }*R^{\\rho}_{\\gamma \\nu \\sigma } +1/20*R^{\\sigma}_{\\beta \\mu \\alpha }*R^{\\rho}_{\\sigma \\nu \\gamma }+1/10*R^{\\sigma}_{\\mu \\beta \\gamma }*R^{\\rho}_{\\nu \\alpha \\sigma } +1/60*R^{\\sigma}_{\\beta \\mu \\gamma }*R^{\\rho}_{\\alpha \\nu \\sigma }+7/120*R_{\\alpha \\beta}*R^{\\rho}_{\\nu \\gamma \\mu } +11/60*R_{\\beta \\mu}*R^{\\rho}_{\\nu \\alpha \\gamma }) *HATK^{\\alpha \\beta  }*DELTA^{\\mu \\nu  }*HATK^{\\gamma  } + L**2*(L-1)*(L-2)*n_{\\rho} *(7/240*R_{\\alpha \\beta}*R^{\\rho}_{\\gamma \\mu \\nu }+7/240*R_{\\alpha \\nu}*R^{\\rho}_{\\beta \\gamma \\mu } -1/60*R_{\\alpha \\mu}*R^{\\rho}_{\\beta \\gamma \\nu }-1/24*R^{\\sigma}_{\\alpha \\beta \\nu }*R^{\\rho}_{\\sigma \\gamma \\mu } +1/15*R^{\\sigma}_{\\alpha \\beta \\nu }*R^{\\rho}_{\\mu \\gamma \\sigma }+1/40*R^{\\sigma}_{\\alpha \\beta \\mu }*R^{\\rho}_{\\sigma \\gamma \\nu } +1/40*R_{\\beta \\gamma}*R^{\\rho}_{\\nu \\mu \\alpha }+1/48*R^{\\sigma}_{\\beta \\gamma \\mu }*R^{\\rho}_{\\nu \\alpha \\sigma }) *HATK^{\\alpha \\beta \\gamma  }*DELTA^{\\mu  }*HATK^{\\nu  } + L**2*(L-1)**2*(L-2) *n_{\\rho}*(-7/240*R^{\\rho}_{\\beta \\gamma \\nu }*R_{\\mu \\alpha}+1/240*R^{\\rho}_{\\mu \\alpha \\nu }*R_{\\beta \\gamma} -1/40*R^{\\rho}_{\\nu \\gamma \\sigma }*R^{\\sigma}_{\\alpha \\mu \\beta }) *HATK^{\\mu \\nu  }*HATK^{\\alpha \\beta \\gamma  } + L*(L-1)*(L-2)*(L-3) *(1/180*R_{\\mu \\nu}*R_{\\alpha \\beta}+7/720*R^{\\sigma}_{\\alpha \\beta \\rho }*R^{\\rho}_{\\mu \\nu \\sigma }) *HATK^{\\mu \\nu \\alpha \\beta  }";
    private static final String DELTA_1_ = "DELTA^\\mu=-L*HATK^\\mu";
    private static final String DELTA_2_ = "DELTA^{\\mu\\nu}=-(1/2)*L*(L-1)*HATK^{\\mu\\nu}+Power[L,2]*(1/2)*(HATK^{\\mu }*HATK^{\\nu }+HATK^{\\nu }*HATK^{\\mu })";
    private static final String DELTA_3_ =
            "DELTA^{\\mu\\nu\\alpha}="
                    + "-(1/6)*L*(L-1)*(L-2)*HATK^{\\mu\\nu\\alpha}"
                    + "+(1/2)*Power[L,2]*(L-1)*(1/3)*("
                    + "HATK^{\\mu \\nu }*HATK^{\\alpha }+"
                    + "HATK^{\\alpha \\nu }*HATK^{\\mu }+"
                    + "HATK^{\\mu \\alpha }*HATK^{\\nu })"
                    + "+1/2*Power[L,2]*(L-1)*(1/3)*("
                    + "HATK^{\\alpha }*HATK^{\\mu \\nu }+"
                    + "HATK^{\\mu }*HATK^{\\alpha \\nu }+"
                    + "HATK^{\\nu }*HATK^{\\alpha \\mu })"
                    + "-Power[L,3]*(1/6)*("
                    + "HATK^{\\mu }*HATK^{\\nu }*HATK^{\\alpha }+"
                    + "HATK^{\\mu }*HATK^{\\alpha }*HATK^{\\nu }+"
                    + "HATK^{\\nu }*HATK^{\\alpha }*HATK^{\\mu }+"
                    + "HATK^{\\nu }*HATK^{\\mu }*HATK^{\\alpha }+"
                    + "HATK^{\\alpha }*HATK^{\\mu }*HATK^{\\nu }+"
                    + "HATK^{\\alpha }*HATK^{\\nu }*HATK^{\\mu })";
    private static final String DELTA_4_ =
            "DELTA^{\\mu\\nu\\alpha\\beta}="
                    + "-(1/24)*L*(L-1)*(L-2)*(L-3)*HATK^{\\mu\\nu\\alpha\\beta}"
                    + "+(1/6)*Power[L,2]*(L-1)*(L-2)*(1/4)*("
                    + "HATK^{\\mu \\nu \\alpha }*HATK^{\\beta }+"
                    + "HATK^{\\mu \\nu \\beta }*HATK^{\\alpha }+"
                    + "HATK^{\\beta \\mu \\alpha }*HATK^{\\nu }+"
                    + "HATK^{\\nu \\beta \\alpha }*HATK^{\\mu })"
                    + "+(1/6)*Power[L,2]*(L-1)*(L-2)*(1/4)*("
                    + "HATK^{\\beta }*HATK^{\\mu \\nu \\alpha }+"
                    + "HATK^{\\alpha }*HATK^{\\mu \\nu \\beta }+"
                    + "HATK^{\\mu }*HATK^{\\beta \\nu \\alpha }+"
                    + "HATK^{\\nu }*HATK^{\\beta \\mu \\alpha })"
                    + "+(1/4)*Power[L,2]*Power[(L-1),2]*(1/6)*("
                    + "HATK^{\\mu\\nu}*HATK^{\\alpha\\beta}+"
                    + "HATK^{\\mu\\beta}*HATK^{\\alpha\\nu}+"
                    + "HATK^{\\mu\\alpha}*HATK^{\\nu\\beta}+"
                    + "HATK^{\\alpha\\nu}*HATK^{\\mu\\beta}+"
                    + "HATK^{\\beta\\nu}*HATK^{\\alpha\\mu}+"
                    + "HATK^{\\alpha\\beta}*HATK^{\\mu\\nu})"
                    + "-(1/2)*Power[L,3]*(L-1)*(1/12)*("
                    + "HATK^{\\mu\\nu}*HATK^\\alpha*HATK^\\beta+"
                    + "HATK^{\\mu\\nu}*HATK^\\beta*HATK^\\alpha+"
                    + "HATK^{\\mu\\beta}*HATK^\\alpha*HATK^\\nu+"
                    + "HATK^{\\mu\\beta}*HATK^\\nu*HATK^\\alpha+"
                    + "HATK^{\\mu\\alpha}*HATK^\\nu*HATK^\\beta+"
                    + "HATK^{\\mu\\alpha}*HATK^\\beta*HATK^\\nu+"
                    + "HATK^{\\nu\\alpha}*HATK^\\mu*HATK^\\beta+"
                    + "HATK^{\\nu\\alpha}*HATK^\\beta*HATK^\\mu+"
                    + "HATK^{\\nu\\beta}*HATK^\\alpha*HATK^\\mu+"
                    + "HATK^{\\nu\\beta}*HATK^\\mu*HATK^\\alpha+"
                    + "HATK^{\\alpha\\beta}*HATK^\\mu*HATK^\\nu+"
                    + "HATK^{\\alpha\\beta}*HATK^\\nu*HATK^\\mu)"
                    + "-(1/2)*Power[L,3]*(L-1)*(1/12)*("
                    + "HATK^\\alpha*HATK^{\\mu\\nu}*HATK^\\beta+"
                    + "HATK^\\beta*HATK^{\\mu\\nu}*HATK^\\alpha+"
                    + "HATK^\\alpha*HATK^{\\mu\\beta}*HATK^\\nu+"
                    + "HATK^\\nu*HATK^{\\mu\\beta}*HATK^\\alpha+"
                    + "HATK^\\nu*HATK^{\\mu\\alpha}*HATK^\\beta+"
                    + "HATK^\\beta*HATK^{\\mu\\alpha}*HATK^\\nu+"
                    + "HATK^\\mu*HATK^{\\nu\\alpha}*HATK^\\beta+"
                    + "HATK^\\beta*HATK^{\\nu\\alpha}*HATK^\\mu+"
                    + "HATK^\\alpha*HATK^{\\nu\\beta}*HATK^\\mu+"
                    + "HATK^\\mu*HATK^{\\nu\\beta}*HATK^\\alpha+"
                    + "HATK^\\mu*HATK^{\\alpha\\beta}*HATK^\\nu+"
                    + "HATK^\\nu*HATK^{\\alpha\\beta}*HATK^\\mu)"
                    + "-(1/2)*Power[L,3]*(L-1)*(1/12)*("
                    + "HATK^\\alpha*HATK^\\beta*HATK^{\\mu\\nu}+"
                    + "HATK^\\beta*HATK^\\alpha*HATK^{\\mu\\nu}+"
                    + "HATK^\\alpha*HATK^\\nu*HATK^{\\mu\\beta}+"
                    + "HATK^\\nu*HATK^\\alpha*HATK^{\\mu\\beta}+"
                    + "HATK^\\nu*HATK^\\beta*HATK^{\\mu\\alpha}+"
                    + "HATK^\\beta*HATK^\\nu*HATK^{\\mu\\alpha}+"
                    + "HATK^\\mu*HATK^\\beta*HATK^{\\nu\\alpha}+"
                    + "HATK^\\beta*HATK^\\mu*HATK^{\\nu\\alpha}+"
                    + "HATK^\\alpha*HATK^\\mu*HATK^{\\nu\\beta}+"
                    + "HATK^\\mu*HATK^\\alpha*HATK^{\\nu\\beta}+"
                    + "HATK^\\mu*HATK^\\nu*HATK^{\\alpha\\beta}+"
                    + "HATK^\\nu*HATK^\\mu*HATK^{\\alpha\\beta})"
                    + "+(1/24)*Power[L,4]*("
                    + "HATK^{\\mu}*HATK^{\\nu}*HATK^{\\alpha}*HATK^{\\beta}+"
                    + "HATK^{\\nu}*HATK^{\\mu}*HATK^{\\alpha}*HATK^{\\beta}+"
                    + "HATK^{\\beta}*HATK^{\\nu}*HATK^{\\alpha}*HATK^{\\mu}+"
                    + "HATK^{\\nu}*HATK^{\\beta}*HATK^{\\alpha}*HATK^{\\mu}+"
                    + "HATK^{\\beta}*HATK^{\\mu}*HATK^{\\alpha}*HATK^{\\nu}+"
                    + "HATK^{\\mu}*HATK^{\\beta}*HATK^{\\alpha}*HATK^{\\nu}+"
                    + "HATK^{\\mu}*HATK^{\\nu}*HATK^{\\beta}*HATK^{\\alpha}+"
                    + "HATK^{\\nu}*HATK^{\\mu}*HATK^{\\beta}*HATK^{\\alpha}+"
                    + "HATK^{\\alpha}*HATK^{\\nu}*HATK^{\\beta}*HATK^{\\mu}+"
                    + "HATK^{\\nu}*HATK^{\\alpha}*HATK^{\\beta}*HATK^{\\mu}+"
                    + "HATK^{\\alpha}*HATK^{\\mu}*HATK^{\\beta}*HATK^{\\nu}+"
                    + "HATK^{\\mu}*HATK^{\\alpha}*HATK^{\\beta}*HATK^{\\nu}+"
                    + "HATK^{\\beta}*HATK^{\\nu}*HATK^{\\mu}*HATK^{\\alpha}+"
                    + "HATK^{\\nu}*HATK^{\\beta}*HATK^{\\mu}*HATK^{\\alpha}+"
                    + "HATK^{\\alpha}*HATK^{\\nu}*HATK^{\\mu}*HATK^{\\beta}+"
                    + "HATK^{\\nu}*HATK^{\\alpha}*HATK^{\\mu}*HATK^{\\beta}+"
                    + "HATK^{\\alpha}*HATK^{\\beta}*HATK^{\\mu}*HATK^{\\nu}+"
                    + "HATK^{\\beta}*HATK^{\\alpha}*HATK^{\\mu}*HATK^{\\nu}+"
                    + "HATK^{\\beta}*HATK^{\\mu}*HATK^{\\nu}*HATK^{\\alpha}+"
                    + "HATK^{\\mu}*HATK^{\\beta}*HATK^{\\nu}*HATK^{\\alpha}+"
                    + "HATK^{\\alpha}*HATK^{\\mu}*HATK^{\\nu}*HATK^{\\beta}+"
                    + "HATK^{\\mu}*HATK^{\\alpha}*HATK^{\\nu}*HATK^{\\beta}+"
                    + "HATK^{\\alpha}*HATK^{\\beta}*HATK^{\\nu}*HATK^{\\mu}+"
                    + "HATK^{\\beta}*HATK^{\\alpha}*HATK^{\\nu}*HATK^{\\mu})";
    private static final String ACTION_ = "counterterms = Flat + WR + SR + SSR + FF + FR + RR";
    private final Expression Flat, WR, SR, SSR, FF, FR, RR, DELTA_1, DELTA_2, DELTA_3, DELTA_4, ACTION;

    private OneLoopCounterterms(Expression Flat, Expression WR, Expression SR, Expression SSR, Expression FF, Expression FR, Expression RR, Expression DELTA_1, Expression DELTA_2, Expression DELTA_3, Expression DELTA_4, Expression ACTION) {
        this.Flat = Flat;
        this.WR = WR;
        this.SR = SR;
        this.SSR = SSR;
        this.FF = FF;
        this.FR = FR;
        this.RR = RR;
        this.DELTA_1 = DELTA_1;
        this.DELTA_2 = DELTA_2;
        this.DELTA_3 = DELTA_3;
        this.DELTA_4 = DELTA_4;
        this.ACTION = ACTION;
    }

    /**
     * Returns the Flat counterterms part
     *
     * @return Flat counterterms part
     */
    public Expression Flat() {
        return Flat;
    }

    /**
     * Returns the WR counterterms part
     *
     * @return WR counterterms part
     */
    public Expression WR() {
        return WR;
    }

    /**
     * Returns the SR counterterms part
     *
     * @return SR counterterms part
     */
    public Expression SR() {
        return SR;
    }

    /**
     * Returns the SSR counterterms part
     *
     * @return SSR counterterms part
     */
    public Expression SSR() {
        return SSR;
    }

    /**
     * Returns the FF counterterms part
     *
     * @return FF counterterms part
     */
    public Expression FF() {
        return FF;
    }

    /**
     * Returns the FR counterterms part
     *
     * @return FR counterterms part
     */
    public Expression FR() {
        return FR;
    }

    /**
     * Returns the RR counterterms part
     *
     * @return RR counterterms part
     */
    public Expression RR() {
        return RR;
    }

    /**
     * Return resulting counterterms, i.e. the Flat + WR + SR + SSR + FF + FR + RR.
     * In order to obtain the divergent part of the one loop effective action, one should
     * integrate counterterms over space volume and multiply on 1/(16*\pi^2*(d-4)) factor.
     *
     * @return resulting counterterms
     */
    public Expression getCounterterms() {
        return ACTION;
    }

    /**
     * Returns \Delta^{\mu ...} tensor, where dots mean 'matrix' indices.
     *
     * @return \Delta^{\mu ...} tensor, where dots mean 'matrix' indices.
     */
    public Expression DELTA_1() {
        return DELTA_1;
    }

    /**
     * Returns \Delta^{\mu\nu ...} tensor, where dots mean 'matrix' indices.
     *
     * @return \Delta^{\mu\nu ...} tensor, where dots mean 'matrix' indices.
     */
    public Expression DELTA_2() {
        return DELTA_2;
    }

    /**
     * Returns \Delta^{\mu\nu\alpha ...} tensor, where dots mean 'matrix' indices.
     *
     * @return \Delta^{\mu\nu\alpha ...} tensor, where dots mean 'matrix' indices.
     */
    public Expression DELTA_3() {
        return DELTA_3;
    }

    /**
     * Returns \Delta^{\mu\nu\alpha\beta ...} tensor, where dots mean 'matrix' indices.
     *
     * @return \Delta^{\mu\nu\alpha\beta ...} tensor, where dots mean 'matrix' indices.
     */
    public Expression DELTA_4() {
        return DELTA_4;
    }

    /**
     * This method performs the calculation of the one-loop counterterms.
     * It also prints the interim results to standard output, during the calculation.
     *
     * @param input input parameters container.
     * @return resulting counterterms container.
     */
    public static OneLoopCounterterms calculateOneLoopCounterterms(OneLoopInput input) {
        //Parsing input strings

        //matrices names
        final String[] matrices = new String[]{
                "KINV", "HATK", "HATW", "HATS", "NABLAS",
                "HATN", "HATF", "NABLAF", "HATM", "DELTA",
                "Flat", "FF", "WR", "SR", "SSR", "FR", "RR", "Kn"};

        //F_{\\mu\\nu} type structure
        final StructureOfIndices F_TYPE_STRUCTURE = new StructureOfIndices(IndexType.GreekLower.getType(), 2);
        //matrices indicator for parse preprocessor
        final Indicator<ParseTokenSimpleTensor> matricesIndicator = new Indicator<ParseTokenSimpleTensor>() {

            @Override
            public boolean is(ParseTokenSimpleTensor object) {
                String name = object.name;
                for (String matrix : matrices)
                    if (name.equals(matrix))
                        return true;
                if (name.equals("F") && object.indices.getStructureOfIndices().equals(F_TYPE_STRUCTURE))
                    return true;
                return false;
            }
        };

        int i, matrixIndicesCount = input.getMatrixIndicesCount(), operatorOrder = input.getOperatorOrder();

        //indices to insert
        int upper[] = new int[matrixIndicesCount / 2], lower[] = upper.clone();
        for (i = 0; i < matrixIndicesCount / 2; ++i) {
            upper[i] = IndicesUtils.createIndex(130 + i, IndexType.GreekLower, true);//30 
            lower[i] = IndicesUtils.createIndex(130 + i + matrixIndicesCount / 2, IndexType.GreekLower, false);
        }

        Expression Flat, WR, SR, SSR, FF, FR, RR, DELTA_1, DELTA_2, DELTA_3, DELTA_4, ACTION;

        //preprocessor for Flat, WR, SR, SSR, FF, FR, RR, counterterms
        IndicesInsertion termIndicesInsertion = new IndicesInsertion(
                IndicesFactory.createSimple(null, upper),
                IndicesFactory.createSimple(null, IndicesUtils.getIndicesNames(upper)),
                matricesIndicator);

        Flat = (Expression) Tensors.parse(Flat_, termIndicesInsertion);
        WR = (Expression) Tensors.parse(WR_, termIndicesInsertion);
        SR = (Expression) Tensors.parse(SR_, termIndicesInsertion);
        SSR = (Expression) Tensors.parse(SSR_, termIndicesInsertion);
        FF = (Expression) Tensors.parse(FF_, termIndicesInsertion);
        FR = (Expression) Tensors.parse(FR_, termIndicesInsertion);
        RR = (Expression) Tensors.parse(RR_, termIndicesInsertion);
        ACTION = (Expression) Tensors.parse(ACTION_, termIndicesInsertion);
        Expression[] terms = new Expression[]{Flat, WR, SR, SSR, FF, FR, RR};

        //preprocessor for DELTA_1,2,3,4
        IndicesInsertion deltaIndicesInsertion = new IndicesInsertion(
                IndicesFactory.createSimple(null, upper),
                IndicesFactory.createSimple(null, lower),
                matricesIndicator);

        DELTA_1 = (Expression) Tensors.parse(DELTA_1_, deltaIndicesInsertion);
        DELTA_2 = (Expression) Tensors.parse(DELTA_2_, deltaIndicesInsertion);
        DELTA_3 = (Expression) Tensors.parse(DELTA_3_, deltaIndicesInsertion);
        DELTA_4 = (Expression) Tensors.parse(DELTA_4_, deltaIndicesInsertion);
        Expression[] deltaExpressions = new Expression[]{DELTA_1, DELTA_2, DELTA_3, DELTA_4};

        Expression FSubstitution = input.getF();
        for (Transformation background : input.getRiemannBackground())
            FSubstitution = (Expression) background.transform(FSubstitution);

        //Calculations        
        Expression[] riemansSubstitutions = new Expression[]{
                FSubstitution,
                Tensors.parseExpression("R_{\\mu \\nu}^{\\mu}_{\\alpha} = R_{\\nu\\alpha}"),
                Tensors.parseExpression("R_{\\mu\\nu}^{\\alpha}_{\\alpha}=0"),
                Tensors.parseExpression("F_{\\mu}^{\\mu}^{\\alpha}_{\\beta}=0"),
                Tensors.parseExpression("R_{\\mu\\nu\\alpha\\beta}*R^{\\mu\\alpha\\nu\\beta}=(1/2)*R_{\\mu\\nu\\alpha\\beta}*R^{\\mu\\nu\\alpha\\beta}"),
                Tensors.parseExpression("R_{\\mu\\nu\\alpha\\beta}*R^{\\mu\\nu\\alpha\\beta}=4*R_{\\mu\\nu}*R^{\\mu\\nu}-R*R"),
                Tensors.parseExpression("R_{\\mu}^{\\mu}= R")
        };


        Expression kronecker = (Expression) Tensors.parse("d_{\\mu}^{\\mu}=4");
        Transformation n2 = new SqrSubs(Tensors.parseSimple("n_\\mu")), n2Transformer = new Transformer(TraverseState.Leaving, new Transformation[]{n2});
        Transformation[] common = new Transformation[]{EliminateMetricsTransformation.ELIMINATE_METRICS, n2Transformer, kronecker};
        Transformation[] all = ArraysUtils.addAll(common, riemansSubstitutions);
        Tensor temp;

        //Calculating Delta- tensors
        System.out.println("Evaluating \\Delta- tensors.");

        //DELTA_1,2
        for (i = 0; i < 2; ++i) {
            temp = deltaExpressions[i];
            temp = input.getL().transform(temp);

            for (Expression hatK : input.getHatQuantities(0))
                temp = hatK.transform(temp);
            temp = ExpandTransformation.expand(temp, common);
            for (Transformation tr : common)
                temp = tr.transform(temp);

            deltaExpressions[i] = (Expression) temp;
            System.out.println("delta" + i + " done");
        }
        Tensor[] combinations;
        Expression[] calculatedCombinations;
        //DELTA_3 //todo for particular values of L some combinations can be neglected
        combinations = new Tensor[]{
                Tensors.parse("HATK^{\\mu\\nu\\alpha}", deltaIndicesInsertion),
                Tensors.parse("HATK^{\\mu\\nu}*HATK^{\\alpha}", deltaIndicesInsertion),
                Tensors.parse("HATK^{\\alpha}*HATK^{\\mu\\nu}", deltaIndicesInsertion),
                Tensors.parse("HATK^{\\mu}*HATK^{\\nu}*HATK^{\\alpha}", deltaIndicesInsertion)
        };
        calculatedCombinations = new Expression[combinations.length];
        System.out.println("Delta3:");
        for (i = 0; i < combinations.length; ++i) {
            temp = combinations[i];
//            System.out.println("Delta3: subs" + i);
            for (Expression hatK : input.getHatQuantities(0))
                temp = hatK.transform(temp);
//            System.out.println("Delta3: expand" + i);
            temp = ExpandTransformation.expand(temp, common);
            for (Transformation tr : common)
                temp = tr.transform(temp);
            calculatedCombinations[i] = Tensors.expression(combinations[i], temp);

        }
        temp = DELTA_3;
        temp = input.getL().transform(temp);
        for (Expression t : calculatedCombinations)
            temp = new NaiveSubstitution(t.get(0), t.get(1)).transform(temp);//t.transform(temp);
//        System.out.println("Delta3:expand");
        temp = ExpandTransformation.expand(temp, common);
//        System.out.println("Delta3:subs");
        for (Transformation tr : common)
            temp = tr.transform(temp);
        System.out.println("Delta3:done");
        deltaExpressions[2] = (Expression) temp;
        //DELTA_4  //todo for different L values somecombinations can be neglected
        combinations = new Tensor[]{
                Tensors.parse("HATK^{\\mu\\nu\\alpha\\beta}", deltaIndicesInsertion),
                Tensors.parse("HATK^{\\mu\\nu\\alpha}*HATK^{\\beta}", deltaIndicesInsertion),
                Tensors.parse("HATK^{\\beta}*HATK^{\\mu\\nu\\alpha }", deltaIndicesInsertion),
                Tensors.parse("HATK^{\\alpha\\beta}*HATK^{\\mu\\nu}", deltaIndicesInsertion),
                Tensors.parse("HATK^{\\mu}*HATK^{\\nu}*HATK^{\\alpha\\beta}", deltaIndicesInsertion),
                Tensors.parse("HATK^{\\mu}*HATK^{\\alpha\\beta}*HATK^{\\nu}", deltaIndicesInsertion),
                Tensors.parse("HATK^{\\alpha\\beta}*HATK^{\\mu}*HATK^{\\nu}", deltaIndicesInsertion),
                Tensors.parse("HATK^{\\beta}*HATK^{\\alpha}*HATK^{\\mu}*HATK^{\\nu}", deltaIndicesInsertion)};
        calculatedCombinations = new Expression[combinations.length];
        System.out.println("Delta4:");
        for (i = 0; i < combinations.length; ++i) {
            temp = combinations[i];
//            System.out.println("Delta4: subs " + i);
            for (Expression hatK : input.getHatQuantities(0))
                temp = hatK.transform(temp);
//            System.out.println("Delta4: expand " + i);
            temp = ExpandTransformation.expand(temp, common);
//            System.out.println("Delta4: tr" + i);
            for (Transformation tr : common)
                temp = tr.transform(temp);
            calculatedCombinations[i] = Tensors.expression(combinations[i], temp);
        }
        temp = DELTA_4;
        temp = input.getL().transform(temp);
        for (Expression t : calculatedCombinations)
            temp = new NaiveSubstitution(t.get(0), t.get(1)).transform(temp);//t.transform(temp);
//        System.out.println("Delta4: expand");
        temp = ExpandTransformation.expand(temp, common);
        System.out.println("Delta4: tr");
        for (Transformation tr : common)
            temp = tr.transform(temp);
        deltaExpressions[3] = (Expression) temp;

        System.out.println("Evaluating \\Delta- tensors done. Evaluating action terms.");

        for (i = 0; i < terms.length; ++i) {
            temp = terms[i];
//            System.out.println(temp.get(0));
            temp = input.getL().transform(temp);

            temp = input.getF().transform(temp);
            temp = input.getHatF().transform(temp);
            for (Transformation riemannBackround : input.getRiemannBackground())
                temp = riemannBackround.transform(temp);

            temp = ExpandTransformation.expand(temp, all);//TODO may be redundant
            for (Transformation tr : all)
                temp = tr.transform(temp);

            for (Expression nabla : input.getNablaS())
                temp = nabla.transform(temp);

            temp = input.getF().transform(temp);
            temp = input.getHatF().transform(temp);

            for (Expression kn : input.getKnQuantities())
                temp = kn.transform(temp);
//            System.out.println("kn " + temp.get(0));

            for (Expression[] hatQuantities : input.getHatQuantities())
                for (Expression hatQ : hatQuantities)
                    temp = hatQ.transform(temp);

//            System.out.println("k " + temp.get(0));
            for (Expression delta : deltaExpressions)
                temp = delta.transform(temp);

//            System.out.println("delta " + temp.get(0));
            temp = ExpandTransformation.expand(temp, all);
            for (Transformation tr : all)
                temp = tr.transform(temp);

//            System.out.println("expand " + temp.get(0));

            //todo remove this line after fixing Redberry #42
            temp = ExpandTransformation.expand(temp);

            temp = new Averaging(Tensors.parseSimple("n_\\mu")).transform(temp);

            temp = ExpandTransformation.expand(temp, all);
            for (Transformation tr : all)
                temp = tr.transform(temp);
//            System.out.println("expand " + temp.get(0));

            temp = ExpandTransformation.expand(temp, all);

//            System.out.println("expand " + temp.get(0));

            terms[i] = (Expression) temp;
            System.out.println(temp);
        }


        for (Expression term : terms)
            ACTION = (Expression) term.transform(ACTION);

        System.out.println(ACTION);

        return new OneLoopCounterterms(Flat, WR, SR, SSR, FF, FR, RR, deltaExpressions[0], deltaExpressions[1], deltaExpressions[2], deltaExpressions[3], ACTION);
    }
}
