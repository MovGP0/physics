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
package org.redberry.physics.util;

import java.util.Arrays;
import redberry.core.indexes.IndexesUtils;
import redberry.core.tensor.ContractionStructure;
import redberry.core.tensor.Product;
import redberry.core.tensor.ProductContent;
import redberry.core.tensor.SimpleTensor;
import redberry.core.tensor.Tensor;
import redberry.core.tensor.TensorContraction;
import redberry.core.tensor.TensorIterator;
import redberry.core.tensor.TensorNumber;
import redberry.core.transformation.Transformation;
import redberry.core.utils.IntArrayList;

/**
 *
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
public class SqrSubs implements Transformation {
    private int name, hashCode;

    public SqrSubs(SimpleTensor st) {
        if (st.getIndexes().size() != 1)
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
        short si = content.getStratchIndexByHash(hashCode);
        if (si == -1)
            return tensor;

        TensorContraction contraction = new TensorContraction(si, new long[]{((long) si) << 16});
        short[] sIndexes = content.getStretchIndex(); //For preformance.
        int index = Arrays.binarySearch(sIndexes, si);
        while (index >= 0 && sIndexes[index--] == si);
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
            if (cs.contractions[index].equals(contraction)
                    && ((indexName = st.getIndexes().get(0)) & 0x80000000) == 0)
                list.add(indexName);
        } while (index < sIndexes.length - 1 && sIndexes[++index] == si);
        int[] indexes = list.asArray();
        Arrays.sort(indexes);
        TensorIterator iterator = tensor.iterator();
        Tensor current;
        while (iterator.hasNext()) {
            current = iterator.next();
            if (!(current instanceof SimpleTensor))
                continue;
            SimpleTensor st = (SimpleTensor) current;
            if (st.getName() != name)
                continue;
            if (Arrays.binarySearch(indexes, IndexesUtils.getNameWithType(st.getIndexes().get(0))) >= 0)
                iterator.remove();
        }
        if (product.getElements().isEmpty())
            return TensorNumber.createONE();
        return tensor.equivalent();
    }
}
