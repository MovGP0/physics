package cc.redberry.physics.OneLoopAction;

import cc.redberry.core.tensor.Expression;

public class Test {
    public static void main(String[] args) {
        Expression K_2 = new Expression("K^{\\mu\\nu}_\\alpha^\\beta=g^{\\mu\\nu}*d_{\\alpha}^{\\beta}-\\lambda/2*(g^{\\mu\\beta}*d_\\alpha^\\nu+g^{\\nu\\beta}*d_\\alpha^\\mu)");
        Expression KINV = new Expression("KINV_\\alpha^\\beta=d_\\alpha^\\beta-\\gamma*n_\\alpha*n^\\beta");
        Expression W = new Expression("W^{\\alpha}_{\\beta}=P^{\\alpha}_{\\beta}+\\lambda/2*R^\\alpha_\\beta");
        OneLoop oneLoop = OneLoop.calculate(K_2, KINV, W);
        System.out.println(oneLoop.ACTION);
    }

}
