package core;

import redberry.core.indexes.Indexes;
import redberry.core.context.CC;
import redberry.core.tensor.Tensor;
import redberry.core.tensor.test.TTest;
import redberry.core.utils.TensorUtils;
import static org.junit.Assert.*;

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

    public static void assertParent(Tensor tensor) {
        assertTrue(TensorUtils.testParentConsistent(tensor));
    }

    public static void assertIndexes(Tensor tensor) {
        assertTrue(TensorUtils.testIndexesConsistent(tensor));
    }

    public static void assertIndexesParity(Tensor... tensors) {
        for (int i = 1; i < tensors.length; ++i)
            assertTrue(tensors[0].getIndexes().equalsIgnoreOrder(tensors[i].getIndexes()));
    }

    public static void assertIndexesParity(Indexes... indexeses) {
        for (int i = 1; i < indexeses.length; ++i)
            assertTrue(indexeses[0].equalsIgnoreOrder(indexeses[i]));
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
}
