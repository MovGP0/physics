package cc.redberry.physics.STEP;

import cc.redberry.core.context.CC;
import cc.redberry.core.tensor.Expression;
import cc.redberry.core.tensor.Tensor;

/**
 * Created by IntelliJ IDEA.
 * User: Konstantin_2
 * Date: 07.04.12
 * Time: 20:36
 * To change this template use File | Settings | File Templates.
 */
public class MainTensors {
    public static final String _Flat = "Flat=1/4*HATS*HATS*HATS*HATS-HATW*HATS*HATS+1/2*HATW*HATW+HATS*HATN-HATM+(L-2)*NABLAS_\\mu*HATW^\\mu" +
            "-L*NABLAS_\\mu*HATW*HATK^\\mu+1/3*((L-1)*NABLAS_\\mu^\\mu*HATS*HATS-L*NABLAS_\\mu*HATK^\\mu*HATS*HATS" +
            "-(L-1)*NABLAS_\\mu*HATS*HATS^\\mu+L*NABLAS_\\mu*HATS*HATS*HATK^\\mu)-1/2*NABLAS_\\mu*NABLAS_\\nu*DELTA^{\\mu\\nu}" +
            "-1/4*(L-1)*(L-2)*NABLAS_\\mu*NABLAS_\\nu^{\\mu\\nu}+1/2*L*(L-1)*1/2*(NABLAS_\\mu*NABLAS_{\\nu }" +
            "+NABLAS_{\\nu }*NABLAS_{\\mu })*HATS^\\nu*HATK^\\mu";

    public static final String WR_ = "WR=-1/2*Pow[L,2]*HATW*HATF_{\\mu\\nu}*Kn^\\mu*HATK^\\nu+1/3*L*HATW*HATK^\\alpha*DELTA^{\\mu\\nu}*n_\\sigma*R^\\sigma_{\\mu\\alpha\\nu}" +
            "+1/3*Pow[L,2]*(L-1)*HATW*HATK^{\\mu\\nu}*HATK^\\alpha*n_\\sigma*R^\\sigma_{\\mu\\alpha\\nu}-1/6*(L-2)*(L-3)*HATW^{\\mu\\nu}*R_{\\mu\\nu}";

    public static final String SR_ = "SR=-1/6*Pow[L,2]*(L-1)*HATS*NABLAF_{\\mu\\alpha\\nu}*Kn^{\\mu\\nu}*HATK^\\alpha" +
            "+2/3*L*HATS*NABLAF_{\\mu\\nu\\alpha}*Kn^\\alpha*DELTA^{\\mu\\nu}" +
            "-1/12*(L-1)*(L-2)*(L-3)*HATS^{\\alpha\\mu\\nu}*NABLAR_{\\alpha\\mu\\nu}" +
            "-1/12*Pow[L,2]*(L-1)*(L-2)*HATS*HATK^{\\mu\\nu\\alpha}*HATK^\\beta*n_\\sigma*NABLAR_\\alpha^\\sigma_{\\mu\\beta\\nu}" +
            "+L*(L-1)*HATS*HATK^{\\mu\\nu}*DELTA^{\\alpha\\beta}*n_\\sigma*(5/12*NABLAR_\\alpha^\\sigma_{\\nu\\beta\\mu}" +
            "-1/12*NABLAR^\\sigma_{\\mu\\alpha\\nu\\beta})" +
            "-1/2*L*HATS*HATK^\\beta*DELTA^{\\mu\\nu\\alpha}*n_\\sigma*NABLAR^\\sigma_{\\alpha\\mu\\beta\\nu}";

    public static final String SSR_ = "SSR=-1/2*L*(L-1)*HATS*HATS^\\mu*HATF_{\\mu\\nu}*HATK^{\\nu}+1/2*Pow[L,2]*HATS*HATS*HATF_{\\mu\\nu}*Kn^{\\mu}*HATK^\\nu" +
            "+1/12*(L-1)*(L-2)*HATS*HATS^{\\mu\\nu}*R_{\\mu\\nu}+1/3*L*(L-1)*HATS*HATS^\\mu*HATK^\\nu*R_{\\mu\\nu}" +
            "+1/6*HATS*HATS*DELTA^{\\mu\\nu}*R_{\\mu\\nu}-1/6*L*(L-1)*(L-2)*HATS*HATS^{\\mu\\nu}*HATK^\\alpha*n_\\sigma*R^\\sigma_{\\mu\\alpha\\nu}" +
            "+1/3*(L-1)*HATS*HATS^\\alpha*DELTA^{\\mu\\nu}*n_\\sigma*R^\\sigma_{\\mu\\alpha\\nu}" +
            "-1/3*Pow[L,2]*(L-1)*HATS*HATS*HATK^{\\mu\\nu}*HATK^\\alpha*n_\\sigma*R^\\sigma_{\\mu\\alpha\\nu}" +
            "-1/3*L*HATS*HATS*HATK^\\alpha*DELTA^{\\mu\\nu}*n_\\sigma*R^\\sigma_{\\mu\\alpha\\nu}";

    public static final String FF_ = "FF=-1/24*Pow[L,2]*Pow[(L-1),2]*HATK^{\\mu\\nu}*F_{\\mu\\alpha}*HATK^{\\alpha\\beta}*F_{\\nu\\beta}" +
            "+1/24*Pow[L,2]*HATK^\\mu*F_{\\beta\\nu}*DELTA^{\\alpha\\beta}*HATK^\\nu*F_{\\alpha\\mu}" +
            "-5/24*Pow[L,2]*HATK^\\mu*F_{\\beta\\mu}*DELTA^{\\alpha\\beta}*HATK^\\nu*F_{\\alpha\\nu}" +
            "-1/48*Pow[L,2]*(L-1)*HATK^\\mu*F_{\\beta\\nu}*DELTA^\\nu*HATK^{\\alpha\\beta}*F_{\\alpha\\mu}" +
            "-1/48*Pow[L,2]*(L-1)*HATK^\\mu*F_{\\beta\\mu}*DELTA^\\nu*HATK^{\\alpha\\beta}*F_{\\alpha\\nu}";

    public static final String FR_ = "FR=1/40*Pow[L,2]*(L-1)*(L-2)*DELTA^\\mu*HATK^\\nu*HATK^{\\alpha\\beta\\gamma}*F_{\\mu\\alpha}*n_\\sigma*R^\\sigma_{\\gamma\\beta\\nu}" +
            "-Pow[L,2]*(L-1)*(L-2)*DELTA^\\nu*HATK^{\\alpha\\beta\\gamma}*HATK^\\mu*n_\\sigma*(1/60*R^\\sigma_{\\beta\\gamma\\mu}*F_{\\alpha\\nu}" +
            "+1/12*R^\\sigma_{\\beta\\gamma\\nu}*F_{\\alpha\\mu})" +
            "+Pow[L,2]*Pow[(L-1),2]*DELTA^\\alpha*HATK^{\\beta\\gamma}*HATK^{\\mu\\nu}*n_\\sigma*(1/60*R^\\sigma_{\\beta\\mu\\gamma}*F_{\\alpha\\nu}" +
            "+1/20*R^\\sigma_{\\alpha\\mu\\gamma}*F_{\\nu\\beta}+1/15*R^\\sigma_{\\gamma\\mu\\alpha}*F_{\\nu\\beta}" +
            "+1/60*R^\\sigma_{\\mu\\nu\\gamma}*F_{\\alpha\\beta})+Pow[L,2]*(L-1)*DELTA^{\\alpha\\beta}*HATK^{\\gamma\\delta}*HATK^{\\mu}" +
            "*n_\\sigma*(4/15*R^\\sigma_{\\delta\\beta\\gamma}*F_{\\alpha\\mu}-1/30*R^\\sigma_{\\beta\\delta\\alpha}*F_{\\gamma\\mu}" +
            "-1/15*R^\\sigma_{\\alpha\\gamma\\mu}*F_{\\beta\\delta}-1/30*R^\\sigma_{\\gamma\\alpha\\mu}*F_{\\beta\\delta})" +
            "+Pow[L,2]*(L-1)*DELTA^{\\alpha\\beta}*HATK^\\gamma*HATK^{\\mu\\nu}*n_\\sigma*(7/60*R^\\sigma_{\\alpha\\beta\\mu}*F_{\\gamma\\nu}" +
            "-11/60*R^\\sigma_{\\beta\\mu\\gamma}*F_{\\alpha\\nu}+1/5*R^\\sigma_{\\mu\\alpha\\gamma}*F_{\\beta\\nu}" +
            "+1/60*R^\\sigma_{\\mu\\alpha\\nu}*F_{\\gamma\\beta})+Pow[L,2]*DELTA^{\\mu\\alpha\\beta}*HATK^\\gamma*HATK^\\nu*n_\\sigma" +
            "*(7/20*R^\\sigma_{\\alpha\\gamma\\beta}*F_{\\nu\\mu}+1/10*R^\\sigma_{\\alpha\\beta\\nu}*F_{\\gamma\\mu})";

    public static final String RR_ = "RR=1/10*Pow[L,2]*HATK^\\delta*DELTA^{\\mu\\nu\\alpha\\beta}*HATK^\\gamma*n_\\sigma*n_\\rho*R^\\sigma_{\\alpha\\beta\\gamma}*R^\\rho_{\\mu\\nu\\delta}+Pow[L,2]*Pow[(L-1),2]*(L-2)*HATK^{\\beta\\gamma\\delta}*DELTA^\\alpha*HATK^{\\mu\\nu}*n_\\sigma*n_\\rho*(2/45*R^\\rho_{\\alpha\\delta\\nu}*R^\\sigma_{\\beta\\mu\\gamma}-1/120*R^\\rho_{\\delta\\alpha\\nu}*R^\\sigma_{\\beta\\mu\\gamma})+Pow[L,2]*(L-1)*HATK^\\delta*DELTA^{\\alpha\\beta\\gamma}*HATK^{\\mu\\nu}*n_\\sigma*n_\\rho*(-1/10*R^\\rho_{\\mu\\gamma\\nu}*R^\\sigma_{\\alpha\\delta\\beta}+1/15*R^\\rho_{\\delta\\alpha\\nu}*R^\\sigma_{\\beta\\mu\\gamma}+1/60*R^\\rho_{\\beta\\delta\\nu}*R^\\sigma_{\\gamma\\mu\\alpha})+Pow[L,2]*Pow[(L-1),2]*HATK^{\\gamma\\delta}*DELTA^{\\alpha\\beta}*HATK^{\\mu\\nu}*n_\\sigma*n_\\rho*(-1/20*R^\\rho_{\\mu\\beta\\nu}*R^\\sigma_{\\delta\\alpha\\gamma}+1/180*R^\\rho_{\\alpha\\nu\\beta}*R^\\sigma_{\\gamma\\delta\\mu}-7/360*R^\\rho_{\\mu\\gamma\\nu}*R^\\sigma_{\\alpha\\delta\\beta}-1/240*R^\\rho_{\\delta\\beta\\nu}*R^\\sigma_{\\gamma\\alpha\\mu}-1/120*R^\\rho_{\\beta\\gamma\\nu}*R^\\sigma_{\\alpha\\delta\\mu}-1/30*R^\\rho_{\\delta\\beta\\nu}*R^\\sigma_{\\alpha\\gamma\\mu})+Pow[L,2]*(L-1)*(L-2)*HATK^\\delta*DELTA^{\\mu\\nu}*HATK^{\\alpha\\beta\\gamma}*n_\\sigma*n_\\rho*(-1/30*R^\\rho_{\\gamma\\nu\\beta}*R^\\sigma_{\\alpha\\delta\\mu}-1/180*R^\\rho_{\\mu\\gamma\\nu}*R^\\sigma_{\\alpha\\beta\\delta}+1/180*R^\\rho_{\\mu\\gamma\\delta}*R^\\sigma_{\\alpha\\beta\\nu})+Pow[L,2]*Pow[(L-1),2]*(L-2)*HATK^{\\mu\\nu}*DELTA^{\\delta}*HATK^{\\alpha\\beta\\gamma}*n_\\sigma*n_\\rho*(1/45*R^\\rho_{\\mu\\gamma\\nu}*R^\\sigma_{\\alpha\\beta\\delta}-1/80*R^\\rho_{\\beta\\nu\\gamma}*R^\\sigma_{\\mu\\alpha\\delta}+1/90*R^\\rho_{\\beta\\nu\\gamma}*R^\\sigma_{\\delta\\alpha\\mu})+Pow[L,2]*(L-1)*HATK^{\\mu\\nu}*DELTA^{\\alpha\\beta\\gamma}*HATK^\\delta*n_\\sigma*n_\\rho*(7/120*R^\\rho_{\\beta\\gamma\\nu}*R^\\sigma_{\\mu\\alpha\\delta}-3/40*R^\\rho_{\\beta\\gamma\\delta}*R^\\sigma_{\\mu\\alpha\\nu}+1/120*R^\\rho_{\\delta\\gamma\\nu}*R^\\sigma_{\\alpha\\beta\\mu})+Pow[L,2]*(L-1)*(L-2)*HATK^{\\alpha\\beta\\gamma}*DELTA^{\\mu\\nu}*HATK^\\delta*n_\\sigma*n_\\rho*(-1/24*R^\\rho_{\\mu\\gamma\\nu}*R^\\sigma_{\\alpha\\beta\\delta}-1/180*R^\\rho_{\\nu\\gamma\\delta}*R^\\sigma_{\\alpha\\beta\\mu}-1/360*R^\\rho_{\\delta\\gamma\\nu}*R^\\sigma_{\\alpha\\beta\\mu})-1/120*Pow[L,2]*(L-1)*(L-2)*(L-3)*HATK^{\\mu\\nu\\alpha\\beta}*DELTA^{\\delta}*HATK^\\gamma*n_\\sigma*n_\\rho*R^\\rho_{\\alpha\\beta\\gamma}*R^\\sigma_{\\mu\\nu\\delta}-1/80*Pow[L,2]*Pow[(L-1),2]*(L-2)*(L-3)*HATK^{\\alpha\\beta\\gamma\\delta}*HATK^{\\mu\\nu}*n_\\sigma*n_\\rho*R^\\rho_{\\beta\\gamma\\mu}*R^\\sigma_{\\alpha\\delta\\nu}+Pow[L,2]*HATK^\\mu*DELTA^{\\alpha\\beta\\gamma}*HATK^\\nu*n_\\rho*(-1/8*R_{\\beta\\gamma}*R^\\rho_{\\nu\\alpha\\mu}+3/20*R_{\\beta\\gamma}*R^\\rho_{\\mu\\alpha\\nu}+3/40*R_{\\alpha\\mu}*R^\\rho_{\\beta\\gamma\\nu}+1/40*R^\\sigma_{\\beta\\gamma\\mu}*R^\\rho_{\\nu\\alpha\\sigma}-3/20*R^\\sigma_{\\alpha\\beta\\mu}*R^\\rho_{\\gamma\\nu\\sigma}+1/10*R^\\sigma_{\\alpha\\beta\\nu}*R^\\rho_{\\gamma\\mu\\sigma})+Pow[L,2]*(L-1)*HATK^\\gamma*DELTA^{\\alpha\\beta}*HATK^{\\mu\\nu}*n_\\rho*(1/20*R_{\\alpha\\nu}*R^\\rho_{\\gamma\\beta\\mu}+1/20*R_{\\alpha\\gamma}*R^\\rho_{\\mu\\beta\\nu}+1/10*R_{\\alpha\\beta}*R^\\rho_{\\mu\\gamma\\nu}+1/20*R^\\sigma_{\\alpha\\nu\\gamma}*R^\\rho_{\\sigma\\beta\\mu}-1/60*R^\\sigma_{\\mu\\alpha\\nu}*R^\\rho_{\\beta\\sigma\\gamma}+1/10*R^\\sigma_{\\alpha\\beta\\gamma}*R^\\rho_{\\mu\\sigma\\nu}-1/12*R^\\sigma_{\\alpha\\beta\\nu}*R^\\rho_{\\mu\\sigma\\gamma})+Pow[L,2]*Pow[(L-1),2]*HATK^{\\alpha\\beta}*DELTA^{\\gamma}*HATK^{\\mu\\nu}*n_\\rho*(1/60*R_{\\alpha\\mu}*R^\\rho_{\\beta\\nu\\gamma}-1/20*R_{\\alpha\\mu}*R^\\rho_{\\gamma\\nu\\beta}+1/120*R_{\\alpha\\beta}*R^\\rho_{\\mu\\nu\\gamma}+3/40*R_{\\alpha\\gamma}*R^\\rho_{\\nu\\beta\\mu}+1/20*R^\\sigma_{\\gamma\\mu\\alpha}*R^\\rho_{\\nu\\sigma\\beta}+1/120*R^\\sigma_{\\alpha\\mu\\gamma}*R^\\rho_{\\beta\\nu\\sigma}-1/40*R^\\sigma_{\\alpha\\mu\\gamma}*R^\\rho_{\\sigma\\nu\\beta}+1/40*R^\\sigma_{\\alpha\\mu\\beta}*R^\\rho_{\\sigma\\nu\\gamma}-1/20*R^\\sigma_{\\alpha\\mu\\beta}*R^\\rho_{\\gamma\\nu\\sigma}-1/40*R^\\sigma_{\\mu\\beta\\nu}*R^\\rho_{\\gamma\\sigma\\alpha})+Pow[L,2]*(L-1)*HATK^{\\alpha\\beta}*DELTA^{\\mu\\nu}*HATK^{\\gamma}*n_\\rho*(1/20*R^\\sigma_{\\mu\\nu\\beta}*R^\\rho_{\\gamma\\sigma\\alpha}-7/60*R^\\sigma_{\\beta\\mu\\alpha}*R^\\rho_{\\gamma\\nu\\sigma}+1/20*R^\\sigma_{\\beta\\mu\\alpha}*R^\\rho_{\\sigma\\nu\\gamma}+1/10*R^\\sigma_{\\mu\\beta\\gamma}*R^\\rho_{\\nu\\alpha\\sigma}+1/60*R^\\sigma_{\\beta\\mu\\gamma}*R^\\rho_{\\alpha\\nu\\sigma}+7/120*R_{\\alpha\\beta}*R^\\rho_{\\nu\\gamma\\mu}+11/60*R_{\\beta\\mu}*R^\\rho_{\\nu\\alpha\\gamma})+Pow[L,2]*(L-1)*(L-2)*HATK^{\\alpha\\beta\\gamma}*DELTA^{\\mu}*HATK^{\\nu}*n_\\rho*(7/240*R_{\\alpha\\beta}*R^\\rho_{\\gamma\\mu\\nu}+7/240*R_{\\alpha\\nu}*R^\\rho_{\\beta\\gamma\\mu}-1/60*R_{\\alpha\\mu}*R^\\rho_{\\beta\\gamma\\nu}-1/24*R^\\sigma_{\\alpha\\beta\\nu}*R^\\rho_{\\sigma\\gamma\\mu}+1/15*R^\\sigma_{\\alpha\\beta\\nu}*R^\\rho_{\\mu\\gamma\\sigma}+1/40*R^\\sigma_{\\alpha\\beta\\mu}*R^\\rho_{\\sigma\\gamma\\nu}+1/40*R_{\\beta\\gamma}*R^\\rho_{\\nu\\mu\\alpha}+1/48*R^\\sigma_{\\beta\\gamma\\mu}*R^\\rho_{\\nu\\alpha\\sigma})+Pow[L,2]*Pow[(L-1),2]*(L-2)*HATK^{\\alpha\\beta\\gamma}*HATK^{\\mu\\nu}*n_\\rho*(-7/240*R_{\\alpha\\mu}*R^\\rho_{\\beta\\gamma\\nu}+1/240*R_{\\beta\\gamma}*R^\\rho_{\\mu\\alpha\\nu}-1/40*R^\\sigma_{\\alpha\\mu\\beta}*R^\\rho_{\\nu\\gamma\\sigma})+L*(L-1)*(L-2)*(L-3)*HATK^{\\mu\\nu\\alpha\\beta}*(1/180*R_{\\mu\\nu}*R_{\\alpha\\beta}+7/720*R^\\sigma_{\\alpha\\beta\\rho}*R^\\rho_{\\mu\\nu\\sigma})";

    public static final String DELTA_1_ = "DELTA^\\mu=-L*HATK^\\mu";

    public static final String DELTA_2_ = "DELTA^{\\mu\\nu}=-1/2*L*(L-1)*HATK^{\\mu\\nu}+Pow[L,2]*1/2*(HATK^{\\mu }*HATK^{\\nu }+HATK^{\\nu }*HATK^{\\mu })";

    public static final String DELTA_3_ = "DELTA^{\\mu\\nu\\alpha}=-1/6*L*(L-1)*(L-2)*HATK^{\\mu\\nu\\alpha}" +
            "+1/2*Pow[L,2]*(L-1)*1/3*(HATK^{\\mu \\nu }*HATK^{\\alpha }+HATK^{\\alpha \\nu }*HATK^{\\mu }+HATK^{\\mu \\alpha }*HATK^{\\nu })" +
            "+1/2*Pow[L,2]*(L-1)*1/3*(HATK^{\\alpha }*HATK^{\\mu \\nu }+HATK^{\\mu }*HATK^{\\alpha \\nu }+HATK^{\\nu }*HATK^{\\alpha \\mu })" +
            "-Pow[L,3]*1/6*(HATK^{\\mu }*HATK^{\\nu }*HATK^{\\alpha }+HATK^{\\nu }*HATK^{\\alpha }*HATK^{\\mu }+HATK^{\\nu }*HATK^{\\mu }*HATK^{\\alpha }" +
            "+HATK^{\\nu }*HATK^{\\alpha }*HATK^{\\mu }+HATK^{\\alpha }*HATK^{\\mu }*HATK^{\\nu }+HATK^{\\alpha }*HATK^{\\nu }*HATK^{\\mu })";

    public static final String DELTA_4_ = "DELTA^{\\mu\\nu\\alpha\\beta}=-1/24*L*(L-1)*(L-2)*(L-3)*HATK^{\\mu\\nu\\alpha\\beta}+1/6*Pow[L,2]*(L-1)*(L-2)" +
            "*1/4*(HATK^{\\mu \\nu \\alpha }*HATK^{\\beta }+HATK^{\\beta \\mu \\alpha }*HATK^{\\nu }+HATK^{\\nu \\beta \\alpha }*HATK^{\\mu }+HATK^{\\nu \\mu \\alpha }*HATK^{\\beta })" +
            "+1/6*Pow[L,2]*(L-1)*(L-2)" +
            "*1/4*(HATK^{\\beta }*HATK^{\\mu \\nu \\alpha }+HATK^{\\mu }*HATK^{\\beta \\nu \\alpha }+HATK^{\\nu }*HATK^{\\beta \\mu \\alpha })" +
            "+1/4*Pow[L,2]*Pow[(L-1),2]" +
            "*1/6*(HATK^{\\mu\\nu}*HATK^{\\alpha\\beta}+HATK^{\\alpha\\beta}*HATK^{\\mu\\nu}+HATK^{\\alpha\\nu}*HATK^{\\mu\\beta}+HATK^{\\mu\\beta}*HATK^{\\alpha\\nu}+" +
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


    public static final Expression Flat = new Expression(_Flat);
    public static final Expression WR = new Expression(WR_);
    public static final Expression SR = new Expression(SR_);
    public static final Expression SSR = new Expression(SSR_);
    public static final Expression FF = new Expression(FF_);
    public static final Expression FR = new Expression(FR_);
    public static final Expression RR = new Expression(RR_);
    public static final Expression DELTA_1 = new Expression(DELTA_1_);
    public static final Expression DELTA_2 = new Expression(DELTA_2_);
    public static final Expression DELTA_3 = new Expression(DELTA_3_);
    public static final Expression DELTA_4 = new Expression(DELTA_4_);

    public static final Expression HATK_0 =
            new Expression("HATK = KINV*K^{\\mu\\nu}*n_{\\nu}*n_\\mu");
    public static final Expression HATK_1 =
            new Expression("HATK^{\\mu} = KINV*K^{\\mu\\nu}*n_{\\nu}");
    public static final Expression HATK_2 =
            new Expression("HATK^{\\mu\\nu} = KINV*K^{\\mu\\nu}");
    public static final Expression HATK_3 =
            new Expression("HATK^{\\mu\\nu\\alpha} = HATK^{\\mu\\nu\\alpha}");
    public static final Expression HATK_4 =
            new Expression("HATK^{\\mu\\nu\\alpha\\beta} = HATK^{\\mu\\nu\\alpha\\beta}");
    public static final Expression HATF_2 =
            new Expression("HATF^{\\mu\\nu} = KINV*F^{\\mu\\nu}");
    public static final Expression ACTION =
            new Expression("ACTION = Flat + WR + SR + SSR + FF + FR + RR ");
    public static final Expression[] HATKs = new Expression[]{HATK_0, HATK_1, HATK_2, HATK_3, HATK_4, HATF_2};
    public static final Expression[] DELTAs = new Expression[]{DELTA_1, DELTA_2, DELTA_3, DELTA_4};
    public static final Expression[] TERMs = new Expression[]{ACTION, Flat, WR, SR, SSR, FF, FR, RR};
    public static final Expression[] ALL = new Expression[]{Flat, WR, SR, SSR, FF, FR, RR, DELTA_1, DELTA_2, DELTA_3, DELTA_4, HATK_0, HATK_1, HATK_2, HATK_3, HATK_4, HATF_2};

    public static final Tensor[] MATRIX_ACTION = new Tensor[]{

    };

    public static final Tensor[] MATRICES = new Tensor[]{
            CC.parse("KINV"),
            CC.parse("HATK"),
            CC.parse("HATK^{\\mu}"),
            CC.parse("HATK^{\\mu\\nu}"),
            CC.parse("HATK^{\\mu\\nu\\alpha}"),
            CC.parse("HATK^{\\mu\\nu\\alpha\\beta}"),
            CC.parse("HATW"),
            CC.parse("HATW^{\\mu}"),
            CC.parse("HATW^{\\mu\\nu}"),
            CC.parse("HATW^{\\mu\\nu\\alpha}"),
            CC.parse("HATS"),
            CC.parse("HATS^{\\mu}"),
            CC.parse("HATS^{\\mu\\nu}"),
            CC.parse("HATS^{\\mu\\nu\\alpha}"),
            CC.parse("NABLAS"),
            CC.parse("NABLAS^{\\mu}"),
            CC.parse("NABLAS^{\\mu\\nu}"),
            CC.parse("NABLAS^{\\mu\\nu\\alpha}"),
            CC.parse("HATN"),
            CC.parse("HATN^{\\mu}"),
            CC.parse("HATN^{\\mu\\nu}"),
            CC.parse("HATN^{\\mu\\nu\\alpha}"),
            CC.parse("HATF"),
            CC.parse("HATF^{\\mu}"),
            CC.parse("HATF^{\\mu\\nu}"),
            CC.parse("HATF^{\\mu\\nu\\alpha}"),
            CC.parse("NABLAF"),
            CC.parse("NABLAF^{\\mu}"),
            CC.parse("NABLAF^{\\mu\\nu}"),
            CC.parse("NABLAF^{\\mu\\nu\\alpha}"),
            CC.parse("HATM"),
            CC.parse("DELTA"),
            CC.parse("DELTA^{\\mu}"),
            CC.parse("DELTA^{\\mu\\nu}"),
            CC.parse("DELTA^{\\mu\\nu\\alpha}"),
            CC.parse("DELTA^{\\mu\\nu\\alpha\\beta}"),
            CC.parse("Flat"),
            CC.parse("FF"),
            CC.parse("WR"),
            CC.parse("SR"),
            CC.parse("SSR"),
            CC.parse("FR"),
            CC.parse("RR")
    };

}
