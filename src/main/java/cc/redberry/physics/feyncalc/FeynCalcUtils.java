package cc.redberry.physics.feyncalc;

import cc.redberry.core.number.Complex;
import cc.redberry.core.tensor.Expression;
import cc.redberry.core.tensor.SimpleTensor;
import cc.redberry.core.tensor.Tensor;
import cc.redberry.core.tensor.Tensors;

import static cc.redberry.core.tensor.Tensors.*;

/**
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
public final class FeynCalcUtils {

    public static Expression[] setMandelstam(Tensor[][] momentums) {
        checkMandelstamInput(momentums);
        SimpleTensor s = parseSimple("s"), t = parseSimple("t"), u = parseSimple("u");

        Expression[] result = new Expression[10];
        int i;
        // (k1,k1) = m1^2, (k2,k2) = m2^2, (k3,k3) = m3^2, (k4,k4) = m4^2
        for (i = 0; i < 4; ++i)
            result[i] = expression(square((SimpleTensor) momentums[i][0]), pow(momentums[i][1], 2));

        //2(k1, k2) = s - k1^2 - k2^2
        //2(k3, k4) = s - k3^2 - k4^2
        result[i++] = expression(multiply(Complex.TWO, contract((SimpleTensor) momentums[0][0], (SimpleTensor) momentums[1][0])),
                sum(s, negate(sum(pow(momentums[0][1], 2), pow(momentums[1][1], 2)))));
        result[i++] = expression(multiply(Complex.TWO, contract((SimpleTensor) momentums[2][0], (SimpleTensor) momentums[3][0])),
                sum(s, negate(sum(pow(momentums[2][1], 2), pow(momentums[3][1], 2)))));

        //-2(k1, k3) = t - k1^2 - k3^2
        //-2(k2, k4) = t - k2^2 - k4^2
        result[i++] = expression(multiply(Complex.MINUSE_TWO, contract((SimpleTensor) momentums[0][0], (SimpleTensor) momentums[2][0])),
                sum(t, negate(sum(pow(momentums[0][1], 2), pow(momentums[2][1], 2)))));
        result[i++] = expression(multiply(Complex.MINUSE_TWO, contract((SimpleTensor) momentums[1][0], (SimpleTensor) momentums[3][0])),
                sum(t, negate(sum(pow(momentums[1][1], 2), pow(momentums[3][1], 2)))));


        //-2(k1, k4) = u - k1^2 - k4^2
        //-2(k2, k3) = u - k2^2 - k3^2
        result[i++] = expression(multiply(Complex.MINUSE_TWO, contract((SimpleTensor) momentums[0][0], (SimpleTensor) momentums[3][0])),
                sum(u, negate(sum(pow(momentums[0][1], 2), pow(momentums[3][1], 2)))));
        result[i++] = expression(multiply(Complex.MINUSE_TWO, contract((SimpleTensor) momentums[1][0], (SimpleTensor) momentums[2][0])),
                sum(u, negate(sum(pow(momentums[1][1], 2), pow(momentums[2][1], 2)))));

        return result;
    }

    private static void checkMandelstamInput(Tensor[][] momentums) {
        if (momentums.length != 4)
            throw new IllegalArgumentException();
        for (int i = 0; i < 4; ++i) {
            Tensor[] pairs = momentums[i];
            if (pairs.length != 2)
                throw new IllegalArgumentException();
            if (pairs[0].getClass() != SimpleTensor.class || pairs[0].getIndices().size() != 1)
                throw new IllegalArgumentException();
            if (!(pairs[0].getClass() == SimpleTensor.class || pairs[0].getClass() == Complex.class) || pairs[1].getIndices().size() != 0)
                throw new IllegalArgumentException();
            for (int j = 0; j < i; ++j)
                if (((SimpleTensor) momentums[j][0]).getName() == ((SimpleTensor) pairs[0]).getName())
                    throw new IllegalArgumentException();
        }
    }

    private static final Tensor contract(SimpleTensor a, SimpleTensor b) {
        return Tensors.multiply(a, Tensors.simpleTensor(b.getName(), a.getIndices().getInverse()));
    }

    private static final Tensor square(SimpleTensor tensor) {
        return Tensors.multiply(tensor, Tensors.simpleTensor(tensor.getName(), tensor.getIndices().getInverse()));
    }

}
