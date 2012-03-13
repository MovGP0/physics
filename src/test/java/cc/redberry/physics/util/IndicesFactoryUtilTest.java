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
package cc.redberry.physics.util;

import cc.redberry.core.parser.ParserIndices;
import cc.redberry.core.tensor.Tensor;
import org.junit.Ignore;
import org.junit.Test;


/**
 *
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
@Ignore
public class IndicesFactoryUtilTest {
    public IndicesFactoryUtilTest() {
    }

    @Test
    public void testCreate() {
        System.out.println(IndicesFactoryUtil.createIndices(new Tensor[0], ParserIndices.parse("^{\\mu\\nu}_{\\alpha\\beta}")));

    }
    
    @Test
    public void testDump() {
        System.out.println(IndicesFactoryUtil.doubleAndDumpIndices(ParserIndices.parse("^ab_mn")));

    }
}
