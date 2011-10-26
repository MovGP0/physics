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

import redberry.core.context.CC;
import redberry.core.tensor.Tensor;
import redberry.core.tensor.test.TTest;
import redberry.core.transformation.Transformation;
import redberry.core.transformation.substitutions.SubstitutionsFactory;
import redberry.core.utils.Indicator;

/**
 * All hat tensors have prefix "HAT".
 * 
 * @author Stanislav Poslavsky
 */
public class OneLoop {
    public static Tensor COUNTERTERMS = CC.parse("G1");
    public final static Transformation COUNTERTERMS_SUBSTITUTION =
            SubstitutionsFactory.createSubstitution("G1 = "
            + "Flat + WR + SR + SSR + FF + FR + RR ");
    public static Tensor RR = CC.parse("RR");
    //temporary L=2
    public static Tensor RR_SUB = CC.parse("(1/10)*L*L*HATK^{\\delta}*DELTA^{\\mu\\nu\\alpha\\beta}*HATK^{\\gamma}*n_{\\sigma}*n_{\\lambda}*R^{\\sigma}_{\\alpha\\beta\\gamma}*R^{\\lambda}_{\\mu\\nu\\delta} + "
            + "L*L*(L-1)*(L-1)*(L-2)*HATK^{\\beta\\gamma\\delta}*DELTA^{\\alpha}*HATK^{\\mu\\nu}*n_{\\sigma}*n_{\\lambda}*((2/45)*R^{\\lambda}_{\\alpha\\delta\\nu}*R^{\\sigma}_{\\beta\\mu\\gamma}-(1/120)*R^{\\lambda}_{\\delta\\alpha\\nu}*R^{\\sigma}_{\\beta\\mu\\gamma}) + "
            + "L*L*(L-1)*HATK^{\\delta}*DELTA^{\\alpha\\beta\\gamma}*HATK^{\\mu\\nu}*n_{\\sigma}*n_{\\lambda}*(-(1/10)*R^{\\lambda}_{\\mu\\gamma\\nu}*R^{\\sigma}_{\\alpha\\delta\\beta}+(1/15)*R^{\\lambda}_{\\delta\\alpha\\nu}*R^{\\sigma}_{\\beta\\mu\\gamma}+(1/60)*R^{\\lambda}_{\\beta\\delta\\nu}*R^{\\sigma}_{\\gamma\\mu\\alpha})+"
            + "L*L*(L-1)*(L-1)*HATK^{\\gamma\\delta}*DELTA^{\\alpha\\beta}*HATK^{\\mu\\nu}*n_{\\sigma}*n_{\\lambda}*(-(1/20)*R^{\\lambda}_{\\mu\\beta\\nu}*R^{\\sigma}_{\\delta\\alpha\\gamma}+(1/180)*R^{\\lambda}_{\\alpha\\nu\\beta}*R^{\\sigma}_{\\gamma\\delta\\mu}-(7/360)*R^{\\lambda}_{\\mu\\gamma\\nu}*R^{\\sigma}_{\\alpha\\delta\\beta}-(1/240)*R^{\\lambda}_{\\delta\\beta\\nu}*R^{\\sigma}_{\\gamma\\alpha\\mu}-(1/120)*R^{\\lambda}_{\\beta\\gamma\\nu}*R^{\\sigma}_{\\alpha\\delta\\mu}-(1/30)*R^{\\lambda}_{\\delta\\beta\\nu}*R^{\\sigma}_{\\alpha\\gamma\\mu})+"
            + "L*L*(L-1)*HATK^{\\mu\\nu}*DELTA^{\\alpha\\beta\\gamma}*HATK^{\\delta}*n_{\\sigma}*n_{\\lambda}*((7/120)*R^{\\lambda}_{\\beta\\gamma\\nu}*R^{\\sigma}_{\\mu\\alpha\\delta}-(3/40)*R^{\\lambda}_{\\beta\\gamma\\delta}*R^{\\sigma}_{\\mu\\alpha\\nu}+(1/120)*R^{\\lambda}_{\\delta\\gamma\\nu}*R^{\\sigma}_{\\alpha\\beta\\mu})+"
            + "L*L*HATK^{\\mu}*DELTA^{\\alpha\\beta\\gamma}*HATK^{\\nu}*{\\nu}_{\\lambda}*(-(1/8)*R_{\\beta\\gamma}*R^{\\lambda}_{\\nu\\alpha\\mu}+(3/20)*R_{\\beta\\gamma}*R^{\\lambda}_{\\mu\\alpha\\nu}+(3/40)*R_{\\alpha\\mu}*R^{\\lambda}_{\\beta\\gamma\\nu}+(1/40)*R^{\\sigma}_{\\beta\\gamma\\mu}*R^{\\lambda}_{\\nu\\alpha\\sigma}-(3/20)*R^{\\sigma}_{\\alpha\\beta\\mu}*R^{\\lambda}_{\\gamma\\nu\\sigma}+(1/10)*R^{\\sigma}_{\\alpha\\beta\\nu}*R^{\\lambda}_{\\gamma\\mu\\sigma})+"
            + "L*L*(L-1)*HATK^{\\gamma}*DELTA^{\\alpha\\beta}*HATK^{\\mu\\nu}*n_{\\lambda}*((1/20)*R_{\\alpha\\nu}*R^{\\lambda}_{\\gamma\\beta\\mu}+(1/20)*R_{\\alpha\\gamma}*R^{\\lambda}_{\\mu\\beta\\nu}+(1/10)*R_{\\alpha\\beta}*R^{\\lambda}_{\\mu\\gamma\\nu}+(1/20)*R^{\\sigma}_{\\alpha\\nu\\gamma}*R^{\\lambda}_{\\sigma\\beta\\mu}-(1/60)*R^{\\sigma}_{\\mu\\alpha\\nu}*R^{\\lambda}_{\\beta\\sigma\\gamma}+(1/10)*R^{\\sigma}_{\\alpha\\beta\\gamma}*R^{\\lambda}_{\\mu\\sigma\\nu}-(1/12)*R^{\\sigma}_{\\alpha\\beta\\nu}*R^{\\lambda}_{\\mu\\sigma\\gamma})+"
            + "L*L*(L-1)*(L-1)*HATK^{\\alpha\\beta}*DELTA^{\\gamma}*HATK^{\\mu\\nu}*n_{\\lambda}*((1/60)*R_{\\alpha\\mu}*R^{\\lambda}_{\\beta\\nu\\gamma}-(1/20)*R_{\\alpha\\mu}*R^{\\lambda}_{\\gamma\\nu\\beta}+(1/120)*R_{\\alpha\\beta}*R^{\\lambda}_{\\mu\\nu\\gamma}+(3/40)*R_{\\alpha\\gamma}*R^{\\lambda}_{\\nu\\beta\\mu}+(1/20)*R^{\\sigma}_{\\gamma\\mu\\alpha}*R^{\\lambda}_{\\nu\\sigma\\beta}+(1/120)*R^{\\sigma}_{\\alpha\\mu\\gamma}*R^{\\lambda}_{\\beta\\nu\\sigma}-(1/40)*R^{\\sigma}_{\\alpha\\mu\\gamma}*R^{\\lambda}_{\\sigma\\nu\\beta}+(1/40)*R^{\\sigma}_{\\alpha\\mu\\beta}*R^{\\lambda}_{\\sigma\\nu\\gamma}-(1/20)*R^{\\sigma}_{\\alpha\\mu\\beta}*R^{\\lambda}_{\\gamma\\nu\\sigma}-(1/40)*R^{\\sigma}_{\\mu\\beta\\nu}*R^{\\lambda}_{\\gamma\\sigma\\alpha})+"
            + "L*L*(L-1)*HATK^{\\alpha\\beta}*DELTA^{\\mu\\nu}*HATK^{\\gamma}*n_{\\lambda}*((1/20)*R^{\\sigma}_{\\mu\\nu\\beta}*R^{\\lambda}_{\\gamma\\sigma\\alpha}-(7/60)*R^{\\sigma}_{\\beta\\mu\\alpha}*R^{\\lambda}_{\\gamma\\nu\\sigma}+(1/20)*R^{\\sigma}_{\\beta\\mu\\alpha}*R^{\\lambda}_{\\sigma\\nu\\gamma}+(1/10)*R^{\\sigma}_{\\mu\\beta\\gamma}*R^{\\lambda}_{\\nu\\alpha\\sigma}+(1/60)*R^{\\sigma}_{\\mu\\beta\\gamma}*R^{\\lambda}_{\\alpha\\nu\\sigma}+(7/120)*R_{\\alpha\\beta}*R^{\\lambda}_{\\nu\\gamma\\mu}+(11/60)*R_{\\beta\\mu}*R^{\\lambda}_{\\nu\\alpha\\gamma})");
    public final Transformation RR_SUBSTITUTION;
    public static Tensor DELTA_1 = CC.parse("DELTA^{\\mu}");
    //temporary without symmetrization
    public static Tensor DELTA_1_SUB = CC.parse("-L*HATK^{\\mu}");
    public static Tensor DELTA_2 = CC.parse("DELTA^{\\mu\\nu}");
    public static Tensor DELTA_2_SUB = CC.parse("-(1/2)*L*(L-1)*HATK^{\\mu\\nu}+L*L*HATK^{\\mu}*HATK^{\\nu}");
    public static Tensor DELTA_3 = CC.parse("DELTA^{\\mu\\nu\\alpha}");
    public static Tensor DELTA_3_SUB = CC.parse("-(1/6)*L*(L-1)*(L-2)*HATK^{\\mu\\nu\\alpha}+(1/12)*L*L*(L-1)*(HATK^{\\mu\\nu}*HATK^{\\alpha}+HATK^{\\mu\\alpha}*HATK^{\\nu}+HATK^{\\alpha\\nu}*HATK^{\\mu}+HATK^{\\nu\\mu}*HATK^{\\alpha}+HATK^{\\nu\\alpha}*HATK^{\\mu}+HATK^{\\alpha\\mu}*HATK^{\\nu})+(1/2)*L*L*(L-1)*HATK^{\\mu}*HATK^{\\nu}*HATK^{\\alpha}");
    public static Tensor DELTA_4 = CC.parse("DELTA^{\\mu\\nu\\alpha\\beta}");
    //temprory wrong
    public static Tensor DELTA_4_SUB = CC.parse("-(1/6)*L*(L-1)*(L-2)*(L-3)*HATK^{\\mu\\nu\\alpha\\beta}");
    public static Tensor HATK_1 = CC.parse("HATK^{\\mu}");
    public static Tensor HATK_2 = CC.parse("HATK^{\\mu\\nu}");
    public static Tensor HATK_3 = CC.parse("HATK^{\\mu\\nu\\alpha}");
    public static Tensor HATK_4 = CC.parse("HATK^{\\mu\\nu\\alpha\\beta}");
    public static Tensor K_2 = CC.parse("K^{\\mu\\nu}");
    public static Tensor K_2_INV = CC.parse("KINV");
    public static Tensor HATK_1_SUB = CC.parse("KINV*K^{\\mu\\nu}*n_{\\nu}");
    public static Tensor HATK_2_SUB = CC.parse("KINV*K^{\\mu\\nu}");
    public static Tensor HATK_3_SUB = CC.parse("HATK^{\\mu\\nu\\alpha}");
    public static Tensor HATK_4_SUB = CC.parse("HATK^{\\mu\\nu\\alpha\\beta}");
    /*
     * Matrices
     */
    public static Tensor MATRIX_K_2 = CC.parse("K^{ij}_{pq}^{\\mu\\nu}");
    //K^{\\mu\\nu}^{ij}_{pq}
    public static Tensor MATRIX_K_2_SUB = CC.parse("g^{\\mu\\nu}*((1/2)*(d^i_p*d^j_q+d^i_q*d^j_p)-(1/4)*g_{pq}*g^{ij})");
    //KINV^{ij}_{pq}
    public static Tensor MATRIX_K_2_INV = CC.parse("KINV^{ij}_{pq}");
    public static Tensor MATRIX_K_2_INV_SUB = CC.parse("(1/2)*(d^i_p*d^j_q+d^i_q*d^j_p)-(1/4)*g_{pq}*g^{ij}");
    /*
     * Defining substitutions
     */
    public final Transformation DELTA_1_SUBSTITUTION;
    public final Transformation DELTA_2_SUBSTITUTION;
    public final Transformation DELTA_3_SUBSTITUTION;
    public final Transformation DELTA_4_SUBSTITUTION;
    public final Transformation HATK_1_SUBSTITUTION;
    public final Transformation HATK_2_SUBSTITUTION;
    public final Transformation HATK_3_SUBSTITUTION;
    public final Transformation HATK_4_SUBSTITUTION;
    public final Transformation MATRIX_K_2_SUBSTITUTION;
    public final Transformation MATRIX_K_2_INV_SUBSTITUTION;
    public final static Indicator<Tensor> matrices = new Indicator<Tensor>() {
        public boolean is(Tensor tensor) {
            return TTest.testEqualstensorStructure(tensor, K_2)
                    || TTest.testEqualstensorStructure(tensor, K_2_INV)
                    || TTest.testEqualstensorStructure(tensor, HATK_1)
                    || TTest.testEqualstensorStructure(tensor, HATK_2)
                    || TTest.testEqualstensorStructure(tensor, HATK_3)
                    || TTest.testEqualstensorStructure(tensor, HATK_4)
                    || TTest.testEqualstensorStructure(tensor, DELTA_1)
                    || TTest.testEqualstensorStructure(tensor, DELTA_2)
                    || TTest.testEqualstensorStructure(tensor, DELTA_3)
                    || TTest.testEqualstensorStructure(tensor, DELTA_4);
        }
    };

    public OneLoop() {
        this.RR_SUBSTITUTION = SubstitutionsFactory.createSubstitution(RR, RR_SUB);
        this.DELTA_1_SUBSTITUTION = SubstitutionsFactory.createSubstitution(DELTA_1, DELTA_1_SUB);
        this.DELTA_2_SUBSTITUTION = SubstitutionsFactory.createSubstitution(DELTA_2, DELTA_2_SUB);
        this.DELTA_3_SUBSTITUTION = SubstitutionsFactory.createSubstitution(DELTA_3, DELTA_3_SUB);
        this.DELTA_4_SUBSTITUTION = SubstitutionsFactory.createSubstitution(DELTA_4, DELTA_4_SUB);
        this.HATK_1_SUBSTITUTION = SubstitutionsFactory.createSubstitution(HATK_1, HATK_1_SUB);
        this.HATK_2_SUBSTITUTION = SubstitutionsFactory.createSubstitution(HATK_2, HATK_2_SUB);
        this.HATK_3_SUBSTITUTION = SubstitutionsFactory.createSubstitution(HATK_3, HATK_3_SUB);
        this.HATK_4_SUBSTITUTION = SubstitutionsFactory.createSubstitution(HATK_4, HATK_4_SUB);
        this.MATRIX_K_2_SUBSTITUTION = SubstitutionsFactory.createSubstitution(MATRIX_K_2, MATRIX_K_2_SUB);
        this.MATRIX_K_2_INV_SUBSTITUTION = SubstitutionsFactory.createSubstitution(MATRIX_K_2_INV, MATRIX_K_2_INV_SUB);
    }
}
