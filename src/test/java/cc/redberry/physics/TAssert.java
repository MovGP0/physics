package cc.redberry.physics;


import cc.redberry.core.context.CC;
import cc.redberry.core.indices.Indices;
import cc.redberry.core.tensor.Tensor;
import cc.redberry.core.tensor.testing.TTest;
import cc.redberry.core.utils.TensorUtils;

/**
 *
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
public class TAssert {
    public static void assertEquals(Tensor target, Tensor expected) {
        assertTrue(TTest.testEquals(target, expected));
    }

    public static void assertEquals(Tensor target, String expected) {
        assertEquals(target, CC.parse(expected));
    }

    public static void assertParity(Tensor target, Tensor expected) {
        assertTrue(TTest.testParity(target, expected));
    }

    public static void assertParity(Tensor target, String expected) {
        assertParity(target, CC.parse(expected));
    }

    public static void assertOpposite(Tensor target, Tensor expected) {
        assertTrue(TTest.testOpposite(target, expected));
    }

    public static void assertOpposite(Tensor target, String expected) {
        assertOpposite(target, CC.parse(expected));
    }

    public static void assertParent(Tensor tensor) {
        assertTrue(TensorUtils.testParentConsistent(tensor));
    }

    public static void assertIndices(Tensor tensor) {
        assertTrue(TensorUtils.testIndicesConsistent(tensor));
    }

    public static void assertIndices(Tensor[] tensors) {
        for (Tensor tensor : tensors)
            assertTrue(TensorUtils.testIndicesConsistent(tensor));
    }

    public static void assertIndicesParity(Tensor... tensors) {
        for (int i = 1; i < tensors.length; ++i)
            assertTrue(tensors[0].getIndices().equalsIgnoreOrder(tensors[i].getIndices()));
    }

    public static void assertIndicesParity(Indices... indiceses) {
        for (int i = 1; i < indiceses.length; ++i)
            assertTrue(indiceses[0].equalsIgnoreOrder(indiceses[i]));
    }

    public static boolean isEquals(Tensor tensor, String what) {
        return TTest.testEquals(tensor, CC.parse(what));
    }

    public static boolean parity(Tensor tensor, String what) {
        return TTest.testParity(tensor, CC.parse(what));
    }

    public static Tensor _(String tensor) {
        return CC.parse(tensor);
    }

    public static void assertTrue(boolean assertion) {
        org.junit.Assert.assertTrue(assertion);
    }

    public static void assertFalse(boolean assertion) {
        org.junit.Assert.assertFalse(assertion);
    }
}
