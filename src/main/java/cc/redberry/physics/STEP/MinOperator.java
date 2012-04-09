package cc.redberry.physics.STEP;

import cc.redberry.core.tensor.Expression;
import cc.redberry.core.transformations.RenameConflictingIndices;
import cc.redberry.transformation.CalculateNumbers;
import cc.redberry.transformation.Transformer;

/**
 * Created by IntelliJ IDEA.
 * User: Konstantin_2
 * Date: 07.04.12
 * Time: 20:49
 * To change this template use File | Settings | File Templates.
 */
public class MinOperator extends MainTensors {
    public static final Expression Kn = new Expression("Kn=1");

    public static final Expression Kn_1 = new Expression("Kn^\\mu=n^\\mu");
    public static final Expression Kn_2 = new Expression("Kn^{\\mu\\nu}=1/(2*l-1)*(g^{\\mu\\nu}+(2*l-2)*n^\\mu*n^\\nu)");
    public static final Expression Kn_3 =
            new Expression("Kn^{\\mu\\nu\\alpha}=1/(2*l-1)*(g^{\\mu\\nu}*n^\\alpha+g^{\\mu\\alpha}*n^\\nu+g^{\\nu\\alpha}*n^\\mu+(2*l-4)*n^\\mu*n^\\nu*n^\\alpha)");
    public static final Expression Kn_4 =
            new Expression("Kn^{\\mu\\nu\\alpha\\beta}=1/((2*l-1)*(2l-3))*(g^{\\mu\\nu}*g^{\\alpha\\beta}" +
                    "+g^{\\mu\\alpha}*g^{\\nu\\beta}+g^{\\nu\\alpha}*g^{\\mu\\beta}" +
                    "+(2*l-4)*(g^{\\alpha\\beta}*n^\\mu*n^\\nu+g^{\\mu\\nu}*n^\\alpha*n^\\beta+g^{\\alpha\\mu}" +
                    "*n^\\beta*n^\\nu+g^{\\beta\\nu}*n^\\alpha*n^\\mu+g^{\\alpha\\nu}*n^\\beta*n^\\mu" +
                    "+g^{\\beta\\mu}*n^\\alpha*n^\\nu)+(2*l-4)*(2*l-6)*n^\\mu*n^\\nu*n^\\alpha*n^\\beta)");
    public static final Expression K_2 = new Expression("K^{\\mu\\nu}_\\alpha^\\beta=g^{\\mu\\nu}*\\delta_{\\alpha}^{\\beta}-\\lambda/2*(g^{\\mu\\beta}*\\delta_\\alpha^\\nu+g^{\\nu\\beta}*\\delta_\\alpha^\\mu)");
    public static final Expression KINV = new Expression("KINV=1");
    public static final Expression L = new Expression("L = 2");

    public static void main(String[] args) {
        for (Expression ex : ALL)
            ex.eval(new Transformer(RenameConflictingIndices.INSTANCE));
        for (Expression ex : ALL)
            ex.eval(L.asSubstitution(), CalculateNumbers.INSTANCE);

    }
}
