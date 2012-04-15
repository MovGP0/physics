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

import cc.redberry.core.combinatorics.PermutationsGenerator;
import cc.redberry.core.context.CC;
import cc.redberry.core.tensor.Tensor;
import org.junit.Test;

/**
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
public class QEDTest {
    public QEDTest() {
    }

    @Test
    public void test1() {
        QED qed = new QED();
        System.out.println(qed.Action);
    }

        String Flat = "Flat=1/4*Pow[HATS,4]-HATW*Pow[HATS,2]+1/2*Pow[HATW,2]+HATS*HATN-HATM+(L-2)*NABLAS_\\mu*HATW^\\mu" +
                "-L*NABLAS_\\mu*HATW*HATK^\\mu+1/3*((L-1)*NABLAS_\\mu^\\mu*HATS*HATS-L*NABLAS_\\mu*HATK^\\mu*HATS*HATS" +
                "-(L-1)*NABLAS_\\mu*HATS*HATS^\\mu+L*NABLAS_\\mu*HATS*HATS*HATK^\\mu)-1/2*NABLAS_\\mu*NABLAS_\\nu*\\Delta^{\\mu\\nu}" +
                "-1/4*(L-1)*(L-2)*NABLAS_\\mu*NABLAS_\\nu^{\\mu\\nu}+1/2*L*(L-1)*1/2*(NABLAS_\\mu*NABLAS_{\\nu }" +
                "+NABLAS_{\\nu }*NABLAS_{\\mu })*HATS^\\nu*HATK^\\mu";

        String WR = "WR=-1/2*Pow[L,2]*HATW*HATF_{\\mu\\nu}*Kn^\\mu*HATK^\\nu+1/3*L*HATW*HATK^\\alpha*\\Delta^{\\mu\\nu}*n_\\sigma*R^\\sigma_{\\mu\\alpha\\nu}" +
                "+1/3*Pow[L,2]*(L-1)*HATW*HATK^{\\mu\\nu}*HATK^\\alpha*n_\\sigma*R^\\sigma_{\\mu\\alpha\\nu}-1/6*(L-2)*(L-3)*HATW^{\\mu\\nu}*R_{\\mu\\nu}";

        String SR = "SR=-1/6*Pow[L,2]*(L-1)*HATS*NABLAF_{\\mu\\alpha\\nu}*Kn^{\\mu\\nu}*HATK^\\alpha" +
                "+2/3*L*HATS*NABLAF_{\\mu\\nu\\alpha}*Kn^\\alpha*\\Delta^{\\mu\\nu}" +
                "-1/12*(L-1)*(L-2)*(L-3)*HATS^{\\alpha\\mu\\nu}*NABLAR_{\\alpha\\mu\\nu}" +
                "-1/12*Pow[L,2]*(L-1)*(L-2)*HATS*HATK^{\\mu\\nu\\alpha}*HATK^\\beta*n_\\sigma*NABLAR_\\alpha^\\sigma_{\\mu\\beta\\nu}" +
                "+L*(L-1)*HATS*HATK^{\\mu\\nu}*\\Delta^{\\alpha\\beta}*n_\\sigma*(5/12*NABLAR_\\alpha^\\sigma_{\\nu\\beta\\mu}" +
                "-1/12*NABLAR^\\sigma_{\\mu\\alpha\\nu\\beta})" +
                "-1/2*L*HATS*HATK^\\beta*\\Delta^{\\mu\\nu\\alpha}*n_\\sigma*NABLAR^\\sigma_{\\alpha\\mu\\beta\\nu}";

        String SSR = "SSR=-1/2*L*(L-1)*HATS*HATS^\\mu*HATF_{\\mu\\nu}*HATK^{\\nu}+1/2*Pow[L,2]*HATS*HATS*HATF_{\\mu\\nu}*Kn^{\\mu}*HATK^\\nu" +
                "+1/12*(L-1)*(L-2)*HATS*HATS^{\\mu\\nu}*R_{\\mu\\nu}+1/3*L*(L-1)*HATS*HATS^\\mu*HATK^\\nu*R_{\\mu\\nu}" +
                "+1/6*HATS*HATS*\\Delta^{\\mu\\nu}*R_{\\mu\\nu}-1/6*L*(L-1)*(L-2)*HATS*HATS^{\\mu\\nu}*HATK^\\alpha*n_\\sigma*R^\\sigma_{\\mu\\alpha\\nu}" +
                "+1/3*(L-1)*HATS*HATS^\\alpha*\\Delta^{\\mu\\nu}*n_\\sigma*R^\\sigma_{\\mu\\alpha\\nu}" +
                "-1/3*Pow[L,2]*(L-1)*HATS*HATS*HATK^{\\mu\\nu}*HATK^\\alpha*n_\\sigma*R^\\sigma_{\\mu\\alpha\\nu}" +
                "-1/3*L*HATS*HATS*HATK^\\alpha*\\Delta^{\\mu\\nu}*n_\\sigma*R^\\sigma_{\\mu\\alpha\\nu}";

        String FF = "FF=-1/24*Pow[L,2]*Pow[(L-1),2]*HATK^{\\mu\\nu}*F_{\\mu\\alpha}*HATK^{\\alpha\\beta}*F_{\\nu\\beta}" +
                "+1/24*Pow[L,2]*HATK^\\mu*F_{\\beta\\nu}*\\Delta^{\\alpha\\beta}*HATK^\\nu*F_{\\alpha\\mu}" +
                "-5/24*Pow[L,2]*HATK^\\mu*F_{\\beta\\mu}*\\Delta^{\\alpha\\beta}*HATK^\\nu*F_{\\alpha\\nu}" +
                "-1/48*Pow[L,2]*(L-1)*HATK^\\mu*F_{\\beta\\nu}*\\Delta^\\nu*HATK^{\\alpha\\beta}*F_{\\alpha\\mu}" +
                "-1/48*Pow[L,2]*(L-1)*HATK^\\mu*F_{\\beta\\mu}*\\Delta^\\nu*HATK^{\\alpha\\beta}*F_{\\alpha\\nu}";

        String FR = "FR=1/40*Pow[L,2]*(L-1)*(L-2)*\\Delta^\\mu*HATK^\\nu*HATK^{\\alpha\\beta\\gamma}*F_{\\mu\\alpha}*n_\\sigma*R^\\sigma_{\\gamma\\beta\\nu}" +
                "-Pow[L,2]*(L-1)*(L-2)*\\Delta^\\nu*HATK^{\\alpha\\beta\\gamma}*HATK^\\mu*n_\\sigma*(1/60*R^\\sigma_{\\beta\\gamma\\mu}*F_{\\alpha\\nu}" +
                "+1/12*R^\\sigma_{\\beta\\gamma\\nu}*F_{\\alpha\\mu})" +
                "+Pow[L,2]*Pow[(L-1),2]*\\Delta^\\alpha*HATK^{\\beta\\gamma}*HATK^{\\mu\\nu}*n_\\sigma*(1/60*R^\\sigma_{\\beta\\mu\\gamma}*F_{\\alpha\\nu}" +
                "+1/20*R^\\sigma_{\\alpha\\mu\\gamma}*F_{\\nu\\beta}+1/15*R^\\sigma_{\\gamma\\mu\\alpha}*F_{\\nu\\beta}" +
                "+1/60*R^\\sigma_{\\mu\\nu\\gamma}*F_{\\alpha\\beta})+Pow[L,2]*(L-1)*\\Delta^{\\alpha\\beta}*HATK^{\\gamma\\delta}*HATK^{\\mu}" +
                "*n_\\sigma*(4/15*R^\\sigma_{\\delta\\beta\\gamma}*F_{\\alpha\\mu}-1/30*R^\\sigma_{\\beta\\delta\\alpha}*F_{\\gamma\\mu}" +
                "-1/15*R^\\sigma_{\\alpha\\gamma\\mu}*F_{\\beta\\delta}-1/30*R^\\sigma_{\\gamma\\alpha\\mu}*F_{\\beta\\delta})" +
                "+Pow[L,2]*(L-1)*\\Delta^{\\alpha\\beta}*HATK^\\gamma*HATK^{\\mu\\nu}*n_\\sigma*(7/60*R^\\sigma_{\\alpha\\beta\\mu}*F_{\\gamma\\nu}" +
                "-11/60*R^\\sigma_{\\beta\\mu\\gamma}*F_{\\alpha\\nu}+1/5*R^\\sigma_{\\mu\\alpha\\gamma}*F_{\\beta\\nu}" +
                "+1/60*R^\\sigma_{\\mu\\alpha\\nu}*F_{\\gamma\\beta})+Pow[L,2]*\\Delta^{\\mu\\alpha\\beta}*HATK^\\gamma*HATK^\\nu*n_\\sigma" +
                "*(7/20*R^\\sigma_{\\alpha\\gamma\\beta}*F_{\\nu\\mu}+1/10*R^\\sigma_{\\alpha\\beta\\nu}*F_{\\gamma\\mu})";

        String RR = "RR=1/10*Pow[L,2]*HATK^\\delta*\\Delta^{\\mu\\nu\\alpha\\beta}*HATK^\\gamma*n_\\sigma*n_\\rho*R^\\sigma_{\\alpha\\beta\\gamma}*R^\\rho_{\\mu\\nu\\delta}+Pow[L,2]*Pow[(L-1),2]*(L-2)*HATK^{\\beta\\gamma\\delta}*\\Delta^\\alpha*HATK^{\\mu\\nu}*n_\\sigma*n_\\rho*(2/45*R^\\rho_{\\alpha\\delta\\nu}*R^\\sigma_{\\beta\\mu\\gamma}-1/120*R^\\rho_{\\delta\\alpha\\nu}*R^\\sigma_{\\beta\\mu\\gamma})+Pow[L,2]*(L-1)*HATK^\\delta*\\Delta^{\\alpha\\beta\\gamma}*HATK^{\\mu\\nu}*n_\\sigma*n_\\rho*(-1/10*R^\\rho_{\\mu\\gamma\\nu}*R^\\sigma_{\\alpha\\delta\\beta}+1/15*R^\\rho_{\\delta\\alpha\\nu}*R^\\sigma_{\\beta\\mu\\gamma}+1/60*R^\\rho_{\\beta\\delta\\nu}*R^\\sigma_{\\gamma\\mu\\alpha})+Pow[L,2]*Pow[(L-1),2]*HATK^{\\gamma\\delta}*\\Delta^{\\alpha\\beta}*HATK^{\\mu\\nu}*n_\\sigma*n_\\rho*(-1/20*R^\\rho_{\\mu\\beta\\nu}*R^\\sigma_{\\delta\\alpha\\gamma}+1/180*R^\\rho_{\\alpha\\nu\\beta}*R^\\sigma_{\\gamma\\delta\\mu}-7/360*R^\\rho_{\\mu\\gamma\\nu}*R^\\sigma_{\\alpha\\delta\\beta}-1/240*R^\\rho_{\\delta\\beta\\nu}*R^\\sigma_{\\gamma\\alpha\\mu}-1/120*R^\\rho_{\\beta\\gamma\\nu}*R^\\sigma_{\\alpha\\delta\\mu}-1/30*R^\\rho_{\\delta\\beta\\nu}*R^\\sigma_{\\alpha\\gamma\\mu})+Pow[L,2]*(L-1)*(L-2)*HATK^\\delta*\\Delta^{\\mu\\nu}*HATK^{\\alpha\\beta\\gamma}*n_\\sigma*n_\\rho*(-1/30*R^\\rho_{\\gamma\\nu\\beta}*R^\\sigma_{\\alpha\\delta\\mu}-1/180*R^\\rho_{\\mu\\gamma\\nu}*R^\\sigma_{\\alpha\\beta\\delta}+1/180*R^\\rho_{\\mu\\gamma\\delta}*R^\\sigma_{\\alpha\\beta\\nu})+Pow[L,2]*Pow[(L-1),2]*(L-2)*HATK^{\\mu\\nu}*\\Delta^{\\delta}*HATK^{\\alpha\\beta\\gamma}*n_\\sigma*n_\\rho*(1/45*R^\\rho_{\\mu\\gamma\\nu}*R^\\sigma_{\\alpha\\beta\\delta}-1/80*R^\\rho_{\\beta\\nu\\gamma}*R^\\sigma_{\\mu\\alpha\\delta}+1/90*R^\\rho_{\\beta\\nu\\gamma}*R^\\sigma_{\\delta\\alpha\\mu})+Pow[L,2]*(L-1)*HATK^{\\mu\\nu}*\\Delta^{\\alpha\\beta\\gamma}*HATK^\\delta*n_\\sigma*n_\\rho*(7/120*R^\\rho_{\\beta\\gamma\\nu}*R^\\sigma_{\\mu\\alpha\\delta}-3/40*R^\\rho_{\\beta\\gamma\\delta}*R^\\sigma_{\\mu\\alpha\\nu}+1/120*R^\\rho_{\\delta\\gamma\\nu}*R^\\sigma_{\\alpha\\beta\\mu})+Pow[L,2]*(L-1)*(L-2)*HATK^{\\alpha\\beta\\gamma}*\\Delta^{\\mu\\nu}*HATK^\\delta*n_\\sigma*n_\\rho*(-1/24*R^\\rho_{\\mu\\gamma\\nu}*R^\\sigma_{\\alpha\\beta\\delta}-1/180*R^\\rho_{\\nu\\gamma\\delta}*R^\\sigma_{\\alpha\\beta\\mu}-1/360*R^\\rho_{\\delta\\gamma\\nu}*R^\\sigma_{\\alpha\\beta\\mu})-1/120*Pow[L,2]*(L-1)*(L-2)*(L-3)*HATK^{\\mu\\nu\\alpha\\beta}*\\Delta^{\\delta}*HATK^\\gamma*n_\\sigma*n_\\rho*R^\\rho_{\\alpha\\beta\\gamma}*R^\\sigma_{\\mu\\nu\\delta}-1/80*Pow[L,2]*Pow[(L-1),2]*(L-2)*(L-3)*HATK^{\\alpha\\beta\\gamma\\delta}*HATK^{\\mu\\nu}*n_\\sigma*n_\\rho*R^\\rho_{\\beta\\gamma\\mu}*R^\\sigma_{\\alpha\\delta\\nu}+Pow[L,2]*HATK^\\mu*\\Delta^{\\alpha\\beta\\gamma}*HATK^\\nu*n_\\rho*(-1/8*R_{\\beta\\gamma}*R^\\rho_{\\nu\\alpha\\mu}+3/20*R_{\\beta\\gamma}*R^\\rho_{\\mu\\alpha\\nu}+3/40*R_{\\alpha\\mu}*R^\\rho_{\\beta\\gamma\\nu}+1/40*R^\\sigma_{\\beta\\gamma\\mu}*R^\\rho_{\\nu\\alpha\\sigma}-3/20*R^\\sigma_{\\alpha\\beta\\mu}*R^\\rho_{\\gamma\\nu\\sigma}+1/10*R^\\sigma_{\\alpha\\beta\\nu}*R^\\rho_{\\gamma\\mu\\sigma})+Pow[L,2]*(L-1)*HATK^\\gamma*\\Delta^{\\alpha\\beta}*HATK^{\\mu\\nu}*n_\\rho*(1/20*R_{\\alpha\\nu}*R^\\rho_{\\gamma\\beta\\mu}+1/20*R_{\\alpha\\gamma}*R^\\rho_{\\mu\\beta\\nu}+1/10*R_{\\alpha\\beta}*R^\\rho_{\\mu\\gamma\\nu}+1/20*R^\\sigma_{\\alpha\\nu\\gamma}*R^\\rho_{\\sigma\\beta\\mu}-1/60*R^\\sigma_{\\mu\\alpha\\nu}*R^\\rho_{\\beta\\sigma\\gamma}+1/10*R^\\sigma_{\\alpha\\beta\\gamma}*R^\\rho_{\\mu\\sigma\\nu}-1/12*R^\\sigma_{\\alpha\\beta\\nu}*R^\\rho_{\\mu\\sigma\\gamma})+Pow[L,2]*Pow[(L-1),2]*HATK^{\\alpha\\beta}*\\Delta^{\\gamma}*HATK^{\\mu\\nu}*n_\\rho*(1/60*R_{\\alpha\\mu}*R^\\rho_{\\beta\\nu\\gamma}-1/20*R_{\\alpha\\mu}*R^\\rho_{\\gamma\\nu\\beta}+1/120*R_{\\alpha\\beta}*R^\\rho_{\\mu\\nu\\gamma}+3/40*R_{\\alpha\\gamma}*R^\\rho_{\\nu\\beta\\mu}+1/20*R^\\sigma_{\\gamma\\mu\\alpha}*R^\\rho_{\\nu\\sigma\\beta}+1/120*R^\\sigma_{\\alpha\\mu\\gamma}*R^\\rho_{\\beta\\nu\\sigma}-1/40*R^\\sigma_{\\alpha\\mu\\gamma}*R^\\rho_{\\sigma\\nu\\beta}+1/40*R^\\sigma_{\\alpha\\mu\\beta}*R^\\rho_{\\sigma\\nu\\gamma}-1/20*R^\\sigma_{\\alpha\\mu\\beta}*R^\\rho_{\\gamma\\nu\\sigma}-1/40*R^\\sigma_{\\mu\\beta\\nu}*R^\\rho_{\\gamma\\sigma\\alpha})+Pow[L,2]*(L-1)*HATK^{\\alpha\\beta}*\\Delta^{\\mu\\nu}*HATK^{\\gamma}*n_\\rho*(1/20*R^\\sigma_{\\mu\\nu\\beta}*R^\\rho_{\\gamma\\sigma\\alpha}-7/60*R^\\sigma_{\\beta\\mu\\alpha}*R^\\rho_{\\gamma\\nu\\sigma}+1/20*R^\\sigma_{\\beta\\mu\\alpha}*R^\\rho_{\\sigma\\nu\\gamma}+1/10*R^\\sigma_{\\mu\\beta\\gamma}*R^\\rho_{\\nu\\alpha\\sigma}+1/60*R^\\sigma_{\\beta\\mu\\gamma}*R^\\rho_{\\alpha\\nu\\sigma}+7/120*R_{\\alpha\\beta}*R^\\rho_{\\nu\\gamma\\mu}+11/60*R_{\\beta\\mu}*R^\\rho_{\\nu\\alpha\\gamma})+Pow[L,2]*(L-1)*(L-2)*HATK^{\\alpha\\beta\\gamma}*\\Delta^{\\mu}*HATK^{\\nu}*n_\\rho*(7/240*R_{\\alpha\\beta}*R^\\rho_{\\gamma\\mu\\nu}+7/240*R_{\\alpha\\nu}*R^\\rho_{\\beta\\gamma\\mu}-1/60*R_{\\alpha\\mu}*R^\\rho_{\\beta\\gamma\\nu}-1/24*R^\\sigma_{\\alpha\\beta\\nu}*R^\\rho_{\\sigma\\gamma\\mu}+1/15*R^\\sigma_{\\alpha\\beta\\nu}*R^\\rho_{\\mu\\gamma\\sigma}+1/40*R^\\sigma_{\\alpha\\beta\\mu}*R^\\rho_{\\sigma\\gamma\\nu}+1/40*R_{\\beta\\gamma}*R^\\rho_{\\nu\\mu\\alpha}+1/48*R^\\sigma_{\\beta\\gamma\\mu}*R^\\rho_{\\nu\\alpha\\sigma})+Pow[L,2]*Pow[(L-1),2]*(L-2)*HATK^{\\alpha\\beta\\gamma}*HATK^{\\mu\\nu}*n_\\rho*(-7/240*R_{\\alpha\\mu}*R^\\rho_{\\beta\\gamma\\nu}+1/240*R_{\\beta\\gamma}*R^\\rho_{\\mu\\alpha\\nu}-1/40*R^\\sigma_{\\alpha\\mu\\beta}*R^\\rho_{\\nu\\gamma\\sigma})+L*(L-1)*(L-2)*(L-3)*HATK^{\\mu\\nu\\alpha\\beta}*(1/180*R_{\\mu\\nu}*R_{\\alpha\\beta}+7/720*R^\\sigma_{\\alpha\\beta\\rho}*R^\\rho_{\\mu\\nu\\sigma})";

        String DELTA_1 = "-L*HATK^\\mu";

        String DELTA_2 = "-1/2*L*(L-1)*HATK^{\\mu\\nu}+Pow[L,2]*1/2*(HATK^{\\mu }*HATK^{\\nu }+HATK^{\\nu }*HATK^{\\mu })";

        String DELTA_3 = "-1/6*L*(L-1)*(L-2)*HATK^{\\mu\\nu\\alpha}" +
                "+1/2*Pow[L,2]*(L-1)*1/3*(HATK^{\\mu \\nu }*HATK^{\\alpha }+HATK^{\\alpha \\nu }*HATK^{\\mu }+HATK^{\\mu \\alpha }*HATK^{\\nu })" +
                "+1/2*Pow[L,2]*(L-1)*1/3*(HATK^{\\alpha }*HATK^{\\mu \\nu }+HATK^{\\mu }*HATK^{\\alpha \\nu }+HATK^{\\nu }*HATK^{\\alpha \\mu })" +
                "-Pow[L,3]*1/6*(HATK^{\\mu }*HATK^{\\nu }*HATK^{\\alpha }+HATK^{\\nu }*HATK^{\\alpha }*HATK^{\\mu }+HATK^{\\nu }*HATK^{\\mu }*HATK^{\\alpha }" +
                "+HATK^{\\nu }*HATK^{\\alpha }*HATK^{\\mu }+HATK^{\\alpha }*HATK^{\\mu }*HATK^{\\nu }+HATK^{\\alpha }*HATK^{\\nu }*HATK^{\\mu })";

        String DELTA_4 = "-1/24*L*(L-1)*(L-2)*(L-3)*HATK^{\\mu\\nu\\alpha\\beta}+1/6*Pow[L,2]*(L-1)*(L-2)" +
                "*1/4*(HATK^{\\mu \\nu \\alpha }*HATK^{\\beta }+HATK^{\\beta \\mu \\alpha }*HATK^{\\nu }+HATK^{\\nu \\beta \\alpha }*HATK^{\\mu }+HATK^{\\nu \\mu \\alpha }*HATK^{\\beta })" +
                "+1/6*Pow[L,2]*(L-1)*(L-2)" +
                "*1/4*(HATK^{\\beta }*HATK^{\\mu \\nu \\alpha }+HATK^{\\mu }*HATK^{\\beta \\nu \\alpha }+HATK^{\\nu }*HATK^{\\beta \\mu \\alpha })" +
                "+1/4*Pow[L,2]*Pow[(L-1),2]" +
                "*1/6*(HATK^{\\mu\\nu}*HATK^{\\alpha\\beta}+HATK^{\\alpha\\beta}*HATK^{\\mu\\nu}+HATK^{\\alpha\\nu}*HATK^{\\mu\\beta}+HATK^{\\mu\\beta}*HATK^{\\alpha\\nu}+\n" +
                "HATK^{\\mu\\alpha}*HATK^{\\nu\\beta}+HATK^{\\beta\\nu}*HATK^{\\alpha\\mu})" +
                "+1/2*Pow[L,3]*(L-1)*" +
                "1/12*(HATK^{\\mu\\nu}*HATK^\\alpha*HATK^\\beta+HATK^{\\alpha\\beta}*HATK^\\mu*HATK^\\nu+HATK^{\\alpha\\nu}*HATK^\\mu*HATK^\\beta+" +
                "HATK^{\\mu\\beta}*HATK^\\alpha*HATK^\\nu+HATK^{\\mu\\alpha}*HATK^\\nu*HATK^\\beta+HATK^{\\beta\\nu}*HATK^\\alpha*HATK^\\mu+" +
                "HATK^{\\mu\\nu}*HATK^\\beta*HATK^\\alpha+HATK^{\\alpha\\beta}*HATK^\\nu*HATK^\\mu+HATK^{\\alpha\\nu}*HATK^\\beta*HATK^\\mu+" +
                "HATK^{\\mu\\beta}*HATK^\\nu*HATK^\\alpha+HATK^{\\mu\\alpha}*HATK^\\beta*HATK^\\nu+HATK^{\\beta\\nu}*HATK^\\mu*HATK^\\alpha)" +
                "+1/2*Pow[L,3]*(L-1)*" +
                "1/12*(HATK^\\alpha*HATK^{\\mu\\nu}*HATK^\\beta+HATK^\\mu*HATK^{\\alpha\\beta}*HATK^\\nu+HATK^\\mu*HATK^{\\alpha\\nu}*HATK^\\beta+" +
                "HATK^\\alpha*HATK^{\\mu\\beta}*HATK^\\nu+HATK^\\nu*HATK^{\\mu\\alpha}*HATK^\\beta+HATK^\\alpha*HATK^{\\beta\\nu}*HATK^\\mu+" +
                "HATK^\\beta*HATK^{\\mu\\nu}*HATK^\\alpha+HATK^\\nu*HATK^{\\alpha\\beta}*HATK^\\mu+HATK^\\beta*HATK^{\\alpha\\nu}*HATK^\\mu+" +
                "HATK^\\nu*HATK^{\\mu\\beta}*HATK^\\alpha+HATK^\\beta*HATK^{\\mu\\alpha}*HATK^\\nu+HATK^\\mu*HATK^{\\beta\\nu}*HATK^\\alpha)" +
                "+1/2*Pow[L,3]*(L-1)*" +
                "1/12*(HATK^\\beta*HATK^\\alpha*HATK^{\\mu\\nu}+HATK^\\nu*HATK^\\mu*HATK^{\\alpha\\beta}+HATK^\\beta*HATK^\\mu*HATK^{\\alpha\\nu}+" +
                "HATK^\\nu*HATK^\\alpha*HATK^{\\mu\\beta}+HATK^\\beta*HATK^\\nu*HATK^{\\mu\\alpha}+HATK^\\mu*HATK^\\alpha*HATK^{\\beta\\nu}+" +
                "HATK^\\alpha*HATK^\\beta*HATK^{\\mu\\nu}+HATK^\\mu*HATK^\\nu*HATK^{\\alpha\\beta}+HATK^\\mu*HATK^\\beta*HATK^{\\alpha\\nu}+" +
                "HATK^\\alpha*HATK^\\nu*HATK^{\\mu\\beta}+HATK^\\nu*HATK^\\beta*HATK^{\\mu\\alpha}+HATK^\\alpha*HATK^\\mu*HATK^{\\beta\\nu})" +
                "+Pow[L,4]*(1/24)*(HATK^\\mu*HATK^\\nu*HATK^\\alpha*HATK^\\beta+HATK^\\mu*HATK^\\nu*HATK^\\beta*HATK^\\alpha+HATK^\\mu*HATK^\\alpha*HATK^\\nu*HATK^\\beta+HATK^\\mu*HATK^\\alpha*HATK^\\beta*HATK^\\nu+" +
                "HATK^\\mu*HATK^\\beta*HATK^\\nu*HATK^\\alpha+HATK^\\mu*HATK^\\beta*HATK^\\alpha*HATK^\\nu+HATK^\\nu*HATK^\\mu*HATK^\\alpha*HATK^\\beta+HATK^\\nu*HATK^\\mu*HATK^\\beta*HATK^\\alpha+" +
                "HATK^\\nu*HATK^\\alpha*HATK^\\mu*HATK^\\beta+HATK^\\nu*HATK^\\alpha*HATK^\\beta*HATK^\\mu+HATK^\\nu*HATK^\\beta*HATK^\\mu*HATK^\\alpha+HATK^\\nu*HATK^\\beta*HATK^\\alpha*HATK^\\mu+" +
                "HATK^\\alpha*HATK^\\mu*HATK^\\nu*HATK^\\beta+HATK^\\alpha*HATK^\\mu*HATK^\\beta*HATK^\\nu+HATK^\\alpha*HATK^\\nu*HATK^\\mu*HATK^\\beta+HATK^\\alpha*HATK^\\nu*HATK^\\beta*HATK^\\mu+" +
                "HATK^\\alpha*HATK^\\beta*HATK^\\mu*HATK^\\nu+HATK^\\alpha*HATK^\\beta*HATK^\\nu*HATK^\\mu+HATK^\\beta*HATK^\\mu*HATK^\\nu*HATK^\\alpha+HATK^\\beta*HATK^\\mu*HATK^\\alpha*HATK^\\nu+" +
                "HATK^\\beta*HATK^\\nu*HATK^\\mu*HATK^\\alpha+HATK^\\beta*HATK^\\nu*HATK^\\alpha*HATK^\\mu+HATK^\\beta*HATK^\\alpha*HATK^\\mu*HATK^\\nu+HATK^\\beta*HATK^\\alpha*HATK^\\nu*HATK^\\mu)";


    String HATW = "\\hat W\\ \\equiv Kn^{-1} W^{\\mu\\nu\\ldots \\alpha} n_\\mu n_\\nu\\ldots n_\\alpha;\n" +
            "\\nonumber\\\\\n" +
            "&&\n" +
            "\\hat K^\\alpha\\ \\equiv Kn^{-1}\n" +
            "K^{\\mu\\nu\\ldots \\beta\\alpha} n_\\mu n_\\nu\\ldots n_\\beta;\n" +
            "\\nonumber\\\\\n" +
            "&&\n" +
            "\\hat S^{\\alpha\\beta}\\ \\equiv Kn^{-1} S^{\\mu\\nu\\ldots\n" +
            "\\gamma\\beta\\alpha} n_\\mu n_\\nu\\ldots n_\\gamma;\\\\\n" +
            "&&\n" +
            "Kn^{-1}{}_i{}^m\\ Kn_m{}^j = \\delta_i{}^j\\nonumber";

    @Test
    public void test2() {
        Tensor _flat = CC.parse(Flat);
        Tensor _WR = CC.parse(WR);
        Tensor _SR = CC.parse(SR);
        Tensor _SSR = CC.parse(SSR);
        Tensor _FF = CC.parse(FR);
        Tensor _FR = CC.parse(FR);
        Tensor _RR = CC.parse(RR);
        Tensor delta_1 = CC.parse(DELTA_1);
        Tensor delat_2 = CC.parse(DELTA_2);
        Tensor delta_3 = CC.parse(DELTA_3);
        Tensor delta_4 = CC.parse(DELTA_4);
        System.out.println(_flat);
    }

    @Test
    public void test3() {
        PermutationsGenerator pg = new PermutationsGenerator(4);
        String[] indices = {"\\mu", "\\nu", "\\alpha", "\\beta"};
        int count = 0;
            while (pg.hasNext()) {
                int[] permutation = pg.next().getPermutation().copy();
                int j = 0;
                for (int i : permutation)
                    System.out.print("HATK^" + indices[i] + (++j != 4? "*" : ""));
                System.out.print("+");
            if (++count %4 == 0)
            System.out.println();
        }

    }

    @Test
    public void HATK_2_2() {
        String[] indices = {"\\mu", "\\nu", "\\alpha", "\\beta"};
        int[][] list = {{0, 1, 2, 3}, {2, 3, 0, 1}, {2, 1, 0, 3}, {0, 3, 2, 1}, {0, 2, 1, 3}, {3, 1, 2, 0}};
        int count = 0;
        String result = "";
        for (int[] l : list) {
            int j = 0;

            for (int i : l) {
                String HATK = "HATK";
                if ((j) % 2 == 0)
                    result += HATK + "^{" + indices[i];
                else
                    result += indices[i] + "}";
                if (j == 1)
                    result += "*";
                ++j;
            }
            result += "+";
            if (++count % 4 == 0)
                System.out.println(result);
        }
    }

    @Test
    public void HATK_1_3() {
        String[] indices = {"\\mu", "\\nu", "\\alpha", "\\beta"};
        int[][] list = {{0, 1, 2, 3}, {2, 3, 0, 1}, {2, 1, 0, 3}, {0, 3, 2, 1}, {0, 2, 1, 3}, {3, 1, 2, 0},
                {0, 1, 3, 2}, {2, 3, 1, 0}, {2, 1, 3, 0}, {0, 3, 1, 2}, {0, 2, 3, 1}, {3, 1, 0, 2}};
        int count = 0;
        String all = "";
        for (int[] l : list) {
            int j = 0;
            String result = "";
            for (int i : l) {
                String HATK = "HATK";
                if (j == 0)
                    result += "*" + HATK + "^{" + indices[i];
                else if (j == 1)
                    result += indices[i] + "}";
                else if (j == 2)
                    result = HATK + "^" + indices[i] + result;
                else if (j == 3)
                    result = HATK + "^" + indices[i] + "*" + result + "+";
                ++j;
            }
            all += result;
            if (++count % 3 == 0) {

                System.out.println(all);
                all = "";
            }
        }
    }
}
