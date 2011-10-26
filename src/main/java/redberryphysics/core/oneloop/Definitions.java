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
import redberry.core.context.ToStringMode;
import redberry.core.tensor.Tensor;
import redberry.core.transformation.Transformation;
import redberry.core.transformation.substitutions.SubstitutionsFactory;

/**
 * All hat tensors have prefix "HAT".
 * 
 * @author Stanislav Poslavsky
 */
public class Definitions {
    public static Tensor COUNTERTERMS = CC.parse("G1");
    public static Tensor RR = CC.parse("RR");
    public final static Transformation COUNTERTERMS_SUBSTITUTION =
            SubstitutionsFactory.createSubstitution("G1 = "
            + "Flat + WR + SR + SSR + FF + FR + RR ");
    //temporary L=2
    public final static Transformation RR_SUBSTITUTION =
            SubstitutionsFactory.createSubstitution("RR = "
            + "(1/10)*L*L*HATK^{\\delta}*DELTA^{\\mu\\nu\\alpha\\beta}*HATK^{\\gamma}*n_{\\sigma}*n_{\\lambda}*R^{\\sigma}_{\\alpha\\beta\\gamma}*R^{\\lambda}_{\\mu\\nu\\delta} + "
            + "L*L*(L-1)*(L-1)*(L-2)*HATK^{\\beta\\gamma\\delta}*DELTA^{\\alpha}*HATK^{\\mu\\nu}*n_{\\sigma}*n_{\\lambda}*((2/45)*R^{\\lambda}_{\\alpha\\delta\\nu}*R^{\\sigma}_{\\beta\\mu\\gamma}-(1/120)*R^{\\lambda}_{\\delta\\alpha\\nu}*R^{\\sigma}_{\\beta\\mu\\gamma}) + "
            + "L*L*(L-1)*HATK^{\\delta}*DELTA^{\\alpha\\beta\\gamma}*HATK^{\\mu\\nu}*n_{\\sigma}*n_{\\lambda}*(-(1/10)*R^{\\lambda}_{\\mu\\gamma\\nu}*R^{\\sigma}_{\\alpha\\delta\\beta}+(1/15)*R^{\\lambda}_{\\delta\\alpha\\nu}*R^{\\sigma}_{\\beta\\mu\\gamma}+(1/60)*R^{\\lambda}_{\\beta\\delta\\nu}*R^{\\sigma}_{\\gamma\\mu\\alpha})+"
            + "L*L*(L-1)*(L-1)*HATK^{\\gamma\\delta}*DELTA^{\\alpha\\beta}*HATK^{\\mu\\nu}*n_{\\sigma}*n_{\\lambda}*(-(1/20)*R^{\\lambda}_{\\mu\\beta\\nu}*R^{\\sigma}_{\\delta\\alpha\\gamma}+(1/180)*R^{\\lambda}_{\\alpha\\nu\\beta}*R^{\\sigma}_{\\gamma\\delta\\mu}-(7/360)*R^{\\lambda}_{\\mu\\gamma\\nu}*R^{\\sigma}_{\\alpha\\delta\\beta}-(1/240)*R^{\\lambda}_{\\delta\\beta\\nu}*R^{\\sigma}_{\\gamma\\alpha\\mu}-(1/120)*R^{\\lambda}_{\\beta\\gamma\\nu}*R^{\\sigma}_{\\alpha\\delta\\mu}-(1/30)*R^{\\lambda}_{\\delta\\beta\\nu}*R^{\\sigma}_{\\alpha\\gamma\\mu})+"
            + "L*L*(L-1)*HATK^{\\mu\\nu}*DELTA^{\\alpha\\beta\\gamma}*HATK^{\\delta}*n_{\\sigma}*n_{\\lambda}*((7/120)*R^{\\lambda}_{\\beta\\gamma\\nu}*R^{\\sigma}_{\\mu\\alpha\\delta}-(3/40)*R^{\\lambda}_{\\beta\\gamma\\delta}*R^{\\sigma}_{\\mu\\alpha\\nu}+(1/120)*R^{\\lambda}_{\\delta\\gamma\\nu}*R^{\\sigma}_{\\alpha\\beta\\mu})+"
            + "L*L*HATK^{\\mu}*DELTA^{\\alpha\\beta\\gamma}*HATK^{\\nu}*{\\nu}_{\\lambda}*(-(1/8)*R_{\\beta\\gamma}*R^{\\lambda}_{\\nu\\alpha\\mu}+(3/20)*R_{\\beta\\gamma}*R^{\\lambda}_{\\mu\\alpha\\nu}+(3/40)*R_{\\alpha\\mu}*R^{\\lambda}_{\\beta\\gamma\\nu}+(1/40)*R^{\\sigma}_{\\beta\\gamma\\mu}*R^{\\lambda}_{\\nu\\alpha\\sigma}-(3/20)*R^{\\sigma}_{\\alpha\\beta\\mu}*R^{\\lambda}_{\\gamma\\nu\\sigma}+(1/10)*R^{\\sigma}_{\\alpha\\beta\\nu}*R^{\\lambda}_{\\gamma\\mu\\sigma})+"
            + "L*L*(L-1)*HATK^{\\gamma}*DELTA^{\\alpha\\beta}*HATK^{\\mu\\nu}*n_{\\lambda}*((1/20)*R_{\\alpha\\nu}*R^{\\lambda}_{\\gamma\\beta\\mu}+(1/20)*R_{\\alpha\\gamma}*R^{\\lambda}_{\\mu\\beta\\nu}+(1/10)*R_{\\alpha\\beta}*R^{\\lambda}_{\\mu\\gamma\\nu}+(1/20)*R^{\\sigma}_{\\alpha\\nu\\gamma}*R^{\\lambda}_{\\sigma\\beta\\mu}-(1/60)*R^{\\sigma}_{\\mu\\alpha\\nu}*R^{\\lambda}_{\\beta\\sigma\\gamma}+(1/10)*R^{\\sigma}_{\\alpha\\beta\\gamma}*R^{\\lambda}_{\\mu\\sigma\\nu}-(1/12)*R^{\\sigma}_{\\alpha\\beta\\nu}*R^{\\lambda}_{\\mu\\sigma\\gamma})+"
            + "L*L*(L-1)*(L-1)*HATK^{\\alpha\\beta}*DELTA^{\\gamma}*HATK^{\\mu\\nu}*n_{\\lambda}*((1/60)*R_{\\alpha\\mu}*R^{\\lambda}_{\\beta\\nu\\gamma}-(1/20)*R_{\\alpha\\mu}*R^{\\lambda}_{\\gamma\\nu\\beta}+(1/120)*R_{\\alpha\\beta}*R^{\\lambda}_{\\mu\\nu\\gamma}+(3/40)*R_{\\alpha\\gamma}*R^{\\lambda}_{\\nu\\beta\\mu}+(1/20)*R^{\\sigma}_{\\gamma\\mu\\alpha}*R^{\\lambda}_{\\nu\\sigma\\beta}+(1/120)*R^{\\sigma}_{\\alpha\\mu\\gamma}*R^{\\lambda}_{\\beta\\nu\\sigma}-(1/40)*R^{\\sigma}_{\\alpha\\mu\\gamma}*R^{\\lambda}_{\\sigma\\nu\\beta}+(1/40)*R^{\\sigma}_{\\alpha\\mu\\beta}*R^{\\lambda}_{\\sigma\\nu\\gamma}-(1/20)*R^{\\sigma}_{\\alpha\\mu\\beta}*R^{\\lambda}_{\\gamma\\nu\\sigma}-(1/40)*R^{\\sigma}_{\\mu\\beta\\nu}*R^{\\lambda}_{\\gamma\\sigma\\alpha})+"
            + "L*L*(L-1)*HATK^{\\alpha\\beta}*DELTA^{\\mu\\nu}*HATK^{\\gamma}*n_{\\lambda}*((1/20)*R^{\\sigma}_{\\mu\\nu\\beta}*R^{\\lambda}_{\\gamma\\sigma\\alpha}-(7/60)*R^{\\sigma}_{\\beta\\mu\\alpha}*R^{\\lambda}_{\\gamma\\nu\\sigma}+(1/20)*R^{\\sigma}_{\\beta\\mu\\alpha}*R^{\\lambda}_{\\sigma\\nu\\gamma}+(1/10)*R^{\\sigma}_{\\mu\\beta\\gamma}*R^{\\lambda}_{\\nu\\alpha\\sigma}+(1/60)*R^{\\sigma}_{\\mu\\beta\\gamma}*R^{\\lambda}_{\\alpha\\nu\\sigma}+(7/120)*R_{\\alpha\\beta}*R^{\\lambda}_{\\nu\\gamma\\mu}+(11/60)*R_{\\beta\\mu}*R^{\\lambda}_{\\nu\\alpha\\gamma})");
//            + "(1/10)*L*L*HATK^{d}*DELTA^{mnab}*HATK^{g}*n_{s}*n_{r}*R^{s}_{abg}*R^{r}_{mnd} + "
//            + "L*L*(L-1)*(L-1)*(L-2)*HATK^{bgd}*DELTA^{a}*HATK^{mn}*n_{s}*n_{r}*((2/45)*R^{r}_{adn}*R^{s}_{bmg}-(1/120)*R^{r}_{dan}*R^{s}_{bmg}) + "
//            + "L*L*(L-1)*HATK^{d}*DELTA^{abg}*HATK^{mn}*n_{s}*n_{r}*(-(1/10)*R^{r}_{mgn}*R^{s}_{adb}+(1/15)*R^{r}_{dan}*R^{s}_{bmg}+(1/60)*R^{r}_{bdn}*R^{s}_{gma})+"
//            + "L*L*(L-1)*(L-1)*HATK^{gd}*DELTA^{ab}*HATK^{mn}*n_{s}*n_{r}*(-(1/20)*R^{r}_{mbn}*R^{s}_{dag}+(1/180)*R^{r}_{anb}*R^{s}_{gdm}-(7/360)*R^{r}_{mgn}*R^{s}_{adb}-(1/240)*R^{r}_{dbn}*R^{s}_{gam}-(1/120)*R^{r}_{bgn}*R^{s}_{adm}-(1/30)*R^{r}_{dbn}*R^{s}_{agm})+"
//            + "L*L*(L-1)*HATK^{mn}*DELTA^{abg}*HATK^{d}*n_{s}*n_{r}*((7/120)*R^{r}_{bgn}*R^{s}_{mad}-(3/40)*R^{r}_{bgd}*R^{s}_{man}+(1/120)*R^{r}_{dgn}*R^{s}_{abm})+"
//            + "L*L*HATK^{m}*DELTA^{abg}*HATK^{n}*{n}_{r}*(-(1/8)*R_{bg}*R^{r}_{nam}+(3/20)*R_{bg}*R^{r}_{man}+(3/40)*R_{am}*R^{r}_{bgn}+(1/40)*R^{s}_{bgm}*R^{r}_{nas}-(3/20)*R^{s}_{abm}*R^{r}_{gns}+(1/10)*R^{s}_{abn}*R^{r}_{gms})+"
//            + "L*L*(L-1)*HATK^{g}*DELTA^{ab}*HATK^{mn}*n_{r}*((1/20)*R_{an}*R^{r}_{gbm}+(1/20)*R_{ag}*R^{r}_{mbn}+(1/10)*R_{ab}*R^{r}_{mgn}+(1/20)*R^{s}_{ang}*R^{r}_{sbm}-(1/60)*R^{s}_{man}*R^{r}_{bsg}+(1/10)*R^{s}_{abg}*R^{r}_{msn}-(1/12)*R^{s}_{abn}*R^{r}_{msg})+"
//            + "L*L*(L-1)*(L-1)*HATK^{ab}*DELTA^{g}*HATK^{mn}*n_{r}*((1/60)*R_{am}*R^{r}_{bng}-(1/20)*R_{am}*R^{r}_{gnb}+(1/120)*R_{ab}*R^{r}_{mng}+(3/40)*R_{ag}*R^{r}_{nbm}+(1/20)*R^{s}_{gma}*R^{r}_{nsb}+(1/120)*R^{s}_{amg}*R^{r}_{bns}-(1/40)*R^{s}_{amg}*R^{r}_{snb}+(1/40)*R^{s}_{amb}*R^{r}_{sng}-(1/20)*R^{s}_{amb}*R^{r}_{gns}-(1/40)*R^{s}_{mbn}*R^{r}_{gsa})+"
//            + "L*L*(L-1)*HATK^{ab}*DELTA^{mn}*HATK^{g}*n_{r}*((1/20)*R^{s}_{mnb}*R^{r}_{gsa}-(7/60)*R^{s}_{bma}*R^{r}_{gns}+(1/20)*R^{s}_{bma}*R^{r}_{sng}+(1/10)*R^{s}_{mbg}*R^{r}_{nas}+(1/60)*R^{s}_{mbg}*R^{r}_{ans}+(7/120)*R_{ab}*R^{r}_{ngm}+(11/60)*R_{bm}*R^{r}_{nag})");
    //temporary without symmetrization
    public final static Transformation DELTA_1_SUBSTITUTION =
            SubstitutionsFactory.createSubstitution("DELTA^{\\mu} = "
            + "-L*HATK^{\\mu}");
    public final static Transformation DELTA_2_SUBSTITUTION =
            SubstitutionsFactory.createSubstitution("DELTA^{\\mu\\nu} = "
            + "-(1/2)*L*(L-1)*HATK^{\\mu\\nu}+L*L*HATK^{\\mu}*HATK^{\\nu}");
    public final static Transformation DELTA_3_SUBSTITUTION =
            SubstitutionsFactory.createSubstitution("DELTA^{\\mu\\nu\\alpha} = "
            + "-(1/6)*L*(L-1)*(L-2)*HATK^{\\mu\\nu\\alpha}+(1/12)*L*L*(L-1)*(HATK^{\\mu\\nu}*HATK^{\\alpha}+HATK^{\\mu\\alpha}*HATK^{\\nu}+HATK^{\\alpha\\nu}*HATK^{\\mu}+HATK^{\\nu\\mu}*HATK^{\\alpha}+HATK^{\\nu\\alpha}*HATK^{\\mu}+HATK^{\\alpha\\mu}*HATK^{\\nu})+(1/2)*L*L*(L-1)*HATK^{\\mu}*HATK^{\\nu}*HATK^{\\alpha}");
//            + "-(1/6)*L*(L-1)*(L-2)*HATK^{mna}+(1/12)*L*L*(L-1)*(HATK^{mn}*HATK^{a}+HATK^{ma}*HATK^{n}+HATK^{an}*HATK^{m}+HATK^{nm}*HATK^{a}+HATK^{na}*HATK^{m}+HATK^{am}*HATK^{n})+(1/2)*L*L*(L-1)*HATK^m*HATK^n*HATK^a");
    //temprory wrong
    public final static Transformation DELTA_4_SUBSTITUTION =
            SubstitutionsFactory.createSubstitution("DELTA^{\\mu\\nu\\alpha\\beta} = "
            + "-(1/6)*L*(L-1)*(L-2)*(L-3)*HATK^{\\mu\\nu\\alpha\\beta}");
}
