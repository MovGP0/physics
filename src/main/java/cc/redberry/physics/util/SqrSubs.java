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

import cc.redberry.core.indices.IndicesUtils;
import cc.redberry.core.tensor.*;
import cc.redberry.core.utils.IntArrayList;
import cc.redberry.transformation.Transformation;
import java.util.Arrays;



/**
 *
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
public class SqrSubs implements Transformation {
    private int name, hashCode;

    public SqrSubs(SimpleTensor st) {
        if (st.getIndices().size() != 1)
            throw new IllegalArgumentException();
        name = st.getName();
        hashCode = st.hashCode();
    }

    @Override
    public Tensor transform(Tensor tensor) {
        if (!(tensor instanceof Product))
            return tensor;
        Product product = (Product) tensor;
        ProductContent content = product.getContent();
        ContractionStructure cs = content.getContractionStructure();
        short si = content.getStretchIndexByHash(hashCode);
        if (si == -1)
            return tensor;

        TensorContraction contraction = new TensorContraction(si, new long[]{((long) si) << 16});
        short[] sIndices = content.getStretchIndex(); //For preformance.
        int index = Arrays.binarySearch(sIndices, si);
        while (index >= 0 && sIndices[index--] == si);
        index++;
        IntArrayList list = new IntArrayList();
        do {
            Tensor t = content.get(index);
            if (!(t instanceof SimpleTensor))
                continue;
            SimpleTensor st = (SimpleTensor) t;
            if (st.getName() != name)
                continue;
            int indexName;
            if (cs.get(index).equals(contraction)
                    && ((indexName = st.getIndices().get(0)) & 0x80000000) == 0)
                list.add(indexName);
        } while (index < sIndices.length - 1 && sIndices[++index] == si);
        int[] indices = list.toArray();
        Arrays.sort(indices);
        TensorIterator iterator = tensor.iterator();
        Tensor current;
        while (iterator.hasNext()) {
            current = iterator.next();
            if (!(current instanceof SimpleTensor))
                continue;
            SimpleTensor st = (SimpleTensor) current;
            if (st.getName() != name)
                continue;
            if (Arrays.binarySearch(indices, IndicesUtils.getNameWithType(st.getIndices().get(0))) >= 0)
                iterator.remove();
        }
        if (product.getElements().isEmpty())
            return TensorNumber.createONE();
        return tensor.equivalent();
    }
}