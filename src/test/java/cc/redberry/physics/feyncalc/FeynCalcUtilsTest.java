package cc.redberry.physics.feyncalc;

import cc.redberry.core.tensor.Expression;
import cc.redberry.core.tensor.Tensor;
import org.junit.Test;

import static cc.redberry.core.tensor.Tensors.parse;

/**
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
public class FeynCalcUtilsTest {
    @Test
    public void testSetMandelstam1() {
        Tensor[][] input = new Tensor[][]{
                {parse("k1_i"), parse("m1")}, {parse("k2_i"), parse("m2")}, {parse("k3_i"), parse("m3")}, {parse("k4_i"), parse("m4")}
        };
        Expression[] ma = FeynCalcUtils.setMandelstam(input);
        // (k1,k1) = m1^2
        // (k2,k2) = m2^2
        // (k3,k3) = m3^2
        // (k4,k4) = m4^2
        // 2(k1, k2) = s - k1^2 - k2^2
        // 2(k3, k4) = s - k3^2 - k4^2
        // -2(k1, k3) = t - k1^2 - k3^2
        // -2(k2, k4) = t - k2^2 - k4^2
        // -2(k1, k4) = u - k1^2 - k4^2
        // -2(k2, k3) = u - k2^2 - k3^2

        for(Expression e: ma)
            System.out.println(e);
    }
}
