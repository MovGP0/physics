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
package cc.redberry.physics;

import cc.redberry.concurrent.OutputPortUnsafe;
import cc.redberry.core.context.CC;
import cc.redberry.core.indexmapping.IndexMappingBuffer;
import cc.redberry.core.indexmapping.IndexMappings;
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
        org.junit.Assert.assertTrue(TTest.testEquals(target, expected));
    }

    public static void assertEquals(Tensor target, String expected) {
        assertEquals(target, CC.parse(expected));
    }

    public static void assertParity(Tensor target, Tensor expected) {
        org.junit.Assert.assertTrue(target.toString(), TTest.testParity(target, expected));
    }

    public static void assertParity(Tensor target, String expected) {
        assertParity(target, CC.parse(expected));
    }

    public static void assertParityFalse(Tensor target, Tensor expected) {
        assertFalse(TTest.testParity(target, expected));
    }

    public static void assertParityFalse(Tensor target, String expected) {
        assertParityFalse(target, CC.parse(expected));
    }

    public static void assertParent(Tensor tensor) {
        org.junit.Assert.assertTrue(TensorUtils.testParentConsistent(tensor));
    }

    public static void assertIndices(Tensor tensor) {
        org.junit.Assert.assertTrue(TensorUtils.testIndicesConsistent(tensor));
    }

    public static void assertIndices(Tensor[] tensors) {
        for (Tensor tensor : tensors)
            assertTrue(TensorUtils.testIndicesConsistent(tensor));
    }

    public static void assertIndicesFalse(Tensor tensor) {
        assertFalse(TensorUtils.testIndicesConsistent(tensor));
    }

    public static void assertIndicesParity(Tensor... tensors) {
        for (int i = 1; i < tensors.length; ++i)
            org.junit.Assert.assertTrue(tensors[0].getIndices().equalsIgnoreOrder(tensors[i].getIndices()));
    }

    public static void assertIndicesParity(Indices... indiceses) {
        for (int i = 1; i < indiceses.length; ++i)
            org.junit.Assert.assertTrue(indiceses[0].equalsIgnoreOrder(indiceses[i]));
    }

    public static boolean isEquals(Tensor tensor, String what) {
        return TTest.testEquals(tensor, CC.parse(what));
    }

    public static boolean parity(Tensor tensor, String what) {
        return TTest.testParity(tensor, CC.parse(what));
    }

    public static void assertTrue(boolean condition) {
        org.junit.Assert.assertTrue(condition);
    }

    public static void assertFalse(boolean condition) {
        org.junit.Assert.assertFalse(condition);
    }

    public static void assertEquals(long expected, long actual) {
        org.junit.Assert.assertEquals(expected, actual);
    }

    public static void assertEquals(int expected, int actual) {
        org.junit.Assert.assertEquals(expected, actual);
    }

    public static void assertEquals(String expected, String actual) {
        org.junit.Assert.assertEquals(expected, actual);
    }

    public static Tensor _(String tensor) {
        return CC.parse(tensor);
    }

    public static void soutMappingsOP(Tensor from, Tensor to) {
        final OutputPortUnsafe<IndexMappingBuffer> opu =
                IndexMappings.createPort(from, to, true);
        IndexMappingBuffer buffer;
        int count = 0;
        while ((buffer = opu.take()) != null) {
            System.out.println(buffer);
            count++;
        }
        System.out.println("Total mappings count " + count);
    }

    public static void soutMappingsOP(String from, String to) {
        soutMappingsOP(_(from), _(to));
    }

    public static void assertOpposite(Tensor target, Tensor expected) {
        assertTrue(TTest.testOpposite(target, expected));
    }

    public static void assertOpposite(Tensor target, String expected) {
        assertOpposite(target, CC.parse(expected));
    }
}
