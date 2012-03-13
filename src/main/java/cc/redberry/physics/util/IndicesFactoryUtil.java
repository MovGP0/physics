/*
 *
 * Redberry: symbolic tensor computations library.
 * Copyright (C) 2010-2011  Stanislav Poslavsky <stvlpos@mail.ru>
 *
 * This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
package cc.redberry.physics.util;

import cc.redberry.core.indexgenerator.IndexGenerator;
import cc.redberry.core.indices.Indices;
import cc.redberry.core.indices.IndicesFactory;
import cc.redberry.core.indices.IndicesStructure;
import cc.redberry.core.parser.ParserIndices;
import cc.redberry.core.tensor.Tensor;
import cc.redberry.core.utils.TensorUtils;




/**
 *
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
public class IndicesFactoryUtil {
    public static Indices createIndices(Tensor[] used_expressions, IndicesStructure indicesStructure) {
        IndexGenerator generator = new IndexGenerator();
        int i;
        Indices expIndices = TensorUtils.getAllIndices(used_expressions);
        for (i = 0; i < expIndices.size(); ++i)
            generator.add(expIndices.get(i));

        int[] i_array = new int[indicesStructure.size()];
        byte type;
        for (i = 0; i < indicesStructure.size(); ++i) {
            type = indicesStructure.get(i);
            i_array[i] = (generator.generate((byte) (type & 0x7F))) | (type << 24);
        }
        return IndicesFactory.createSimple(i_array);
    }

    public static Indices createIndices(Tensor[] used_expressions, Indices sample) {
        return createIndices(used_expressions, new IndicesStructure(sample));
    }

    public static Indices createIndices(Tensor[] used_expressions, String sample) {
        return createIndices(used_expressions, ParserIndices.parse(sample));
    }

    public static Indices doubleAndDumpIndices(Indices indices) {
        int length = indices.size();
        int[] i_array = indices.getAllIndices().copy();
        int[] res = new int[length * 2];
        System.arraycopy(i_array, 0, res, 0, length);
        for (int i = 0; i < length; ++i)
            res[i + length] = 0x80000000 ^ i_array[i];
        return IndicesFactory.createSimple(res);
    }

    public static Indices doubleAndDumpIndices(String indices) {
        return doubleAndDumpIndices(ParserIndices.parse(indices));
    }
}
