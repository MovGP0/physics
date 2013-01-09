/*
 * Redberry: symbolic tensor computations.
 *
 * Copyright (c) 2010-2013:
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
package cc.redberry.physics.feyncalc;

import cc.redberry.core.TAssert;
import cc.redberry.core.context.CC;
import cc.redberry.core.indices.IndexType;
import cc.redberry.core.parser.preprocessor.GeneralIndicesInsertion;
import cc.redberry.core.tensor.Tensor;
import cc.redberry.core.transformations.EliminateFromSymmetriesTransformation;
import org.junit.Test;

import static cc.redberry.core.tensor.Tensors.*;

/**
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
public class UnitaryTraceTest {

    @Test
    public void test1() {
        GeneralIndicesInsertion indicesInsertion = new GeneralIndicesInsertion();
        CC.current().getParseManager().defaultParserPreprocessors.add(indicesInsertion);
        indicesInsertion.addInsertionRule(parseSimple("T^a'_b'a"), IndexType.Matrix1);

        Tensor t = parse("Tr[T_a*T_b]");
        t = UnitaryTrace.unitaryTrace(t);
        TAssert.assertEquals(t, parse("g_ab/2"));
    }

    @Test
    public void test2() {
        GeneralIndicesInsertion indicesInsertion = new GeneralIndicesInsertion();
        CC.current().getParseManager().defaultParserPreprocessors.add(indicesInsertion);
        indicesInsertion.addInsertionRule(parseSimple("T^a'_b'a"), IndexType.Matrix1);

        setSymmetric(parseSimple("d_abd"));
        addAntiSymmetry(parseSimple("f_abc"), 1, 0, 2);
        addSymmetry("f_abc", 2, 0, 1);
        Tensor t = parse("Tr[T_a*T_b*T_c]");
        t = UnitaryTrace.unitaryTrace(t);
        t = EliminateFromSymmetriesTransformation.ELIMINATE_FROM_SYMMETRIES.transform(t);
        Tensor expected = parse("d_abc/4+I/4*f_abc");
        TAssert.assertEquals(t, expected);
    }


    @Test
    public void test3() {
        GeneralIndicesInsertion indicesInsertion = new GeneralIndicesInsertion();
        CC.current().getParseManager().defaultParserPreprocessors.add(indicesInsertion);
        indicesInsertion.addInsertionRule(parseSimple("T^a'_b'a"), IndexType.Matrix1);

        setSymmetric(parseSimple("d_abd"));
        addAntiSymmetry(parseSimple("f_abc"), 1, 0, 2);
        addSymmetry("f_abc", 2, 0, 1);
        Tensor t = parse("Tr[T_a*T_b*T_c*T_d]");
        Tensor t1 = parse("Tr[T_b*T_a*T_c*T_d]");

        t = UnitaryTrace.unitaryTrace(t);
        t1 = UnitaryTrace.unitaryTrace(t1);
        System.out.println(subtract(t, t1));

        t = EliminateFromSymmetriesTransformation.ELIMINATE_FROM_SYMMETRIES.transform(t);
        Tensor expected = parse("-(I/8)*f_adx*d_bc^x + (I/8)*d_adx*f_bc^x+1/8*d_ade*d_bc^e - 1/8*d_bde*d_ac^e+1/8*d_cde*d_ab^e + 1/(4*N)*g_ad*g_bc - 1/(4*N)*g_ac*g_bd + 1/(4*N)*g_ab*g_cd");
        System.out.println(t);
        System.out.println(expected);
//        TAssert.assertEquals(t, expected);
    }
}
