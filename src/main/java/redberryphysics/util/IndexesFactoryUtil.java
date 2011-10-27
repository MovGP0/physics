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
package redberryphysics.util;

import redberry.core.indexes.Indexes;
import redberry.core.indexes.IndexesFactory;
import redberry.core.indexes.IndexesStructure;
import redberry.core.indexgenerator.IndexGenerator;
import redberry.core.parser.ParserIndexes;
import redberry.core.tensor.Tensor;

/**
 *
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
public class IndexesFactoryUtil {
    public static Indexes createIndexes(Tensor[] used_expressions, IndexesStructure indexesStructure) {
        IndexGenerator generator = new IndexGenerator();
        Indexes expIndexes;
        int i;
        for (Tensor exp : used_expressions) {
            expIndexes = exp.getIndexes();
            for (i = 0; i < expIndexes.size(); ++i)
                generator.add(expIndexes.get(i));
        }
        int[] i_array = new int[indexesStructure.size()];
        byte type;
        for (i = 0; i < indexesStructure.size(); ++i) {
            type = indexesStructure.get(i);
            i_array[i] = (generator.generate((byte) (type & 0x7F))) | (type << 24);
        }
        return IndexesFactory.create(i_array);
    }

    public static Indexes createIndexes(Tensor[] used_expressions, Indexes sample) {
        return createIndexes(used_expressions, new IndexesStructure(sample));
    }

    public static Indexes createIndexes(Tensor[] used_expressions, String sample) {
        return createIndexes(used_expressions, ParserIndexes.parse(sample));
    }

    public static Indexes doubleAndDumpIndexes(Indexes indexes) {
        int length = indexes.size();
        int[] i_array = indexes.getAllIndexes().copy();
        int[] res = new int[length * 2];
        System.arraycopy(i_array, 0, res, 0, length);
        for (int i = 0; i < length; ++i)
            res[i + length] = 0x80000000 ^ i_array[i];
        return IndexesFactory.create(res);
    }

    public static Indexes doubleAndDumpIndexes(String indexes) {
        return doubleAndDumpIndexes(ParserIndexes.parse(indexes));
    }
}
