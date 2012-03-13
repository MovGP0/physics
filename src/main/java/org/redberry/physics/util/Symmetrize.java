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

import cc.redberry.core.number.ComplexElement;
import cc.redberry.core.number.NumberFraction;
import cc.redberry.core.number.RationalElement;
import cc.redberry.core.tensor.Product;
import cc.redberry.core.tensor.Sum;
import cc.redberry.core.tensor.Tensor;
import cc.redberry.core.tensor.TensorNumber;
import cc.redberry.tensorgenerator.IndexMappingPermutationsGenerator;
import cc.redberry.transformation.Transformation;
import java.util.List;



/**
 *
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
public class Symmetrize implements Transformation {
    public static final Symmetrize INSTANCE = new Symmetrize();

    private Symmetrize() {
    }

    @Override
    public Tensor transform(Tensor tensor) {
        if (tensor instanceof Product)
            return symmetrizeProduct((Product) tensor);
        if (tensor instanceof Sum) {
            Sum res = new Sum();
            for (Tensor t : tensor)
                res.add(symmetrizeProduct((Product) t));
            return res;
        }
        throw new UnsupportedOperationException();
    }

    private static Tensor symmetrizeProduct(Product p) {
        List<Tensor> list = IndexMappingPermutationsGenerator.getAllPermutations(p);
        return new Product(
                new TensorNumber(
                new ComplexElement(
                new NumberFraction((long) 1, (long) list.size()), RationalElement.ZERO)), new Sum(list));
    }
//    private final Indices indices;
//    private final int lowerCount;
//    private final int upperCount;
//    private final int[] indicesNames;
//
//    public Symmetrize(Indices indices) {
//        this.indices = indices.getFreeIndices();
//        lowerCount = indices.getLower().size();
//        upperCount = indices.getUpper().size();
//        indicesNames = new int[this.indices.size()];
//        for (int i = 0; i < indicesNames.length; ++i)
//            indicesNames[i] = IndicesUtils.getNameWithType(this.indices.get(i));
//    }
//
//    public Symmetrize(String indices) {
//        this(ParserIndices.parse(indices));
//    }
//
//    @Override
//    public Tensor transform(Tensor tensor) {
//        if (!tensor.getIndices().getFreeIndices().equalsIgnoreOrder(indices))
//            throw new IllegalArgumentException();
//
//        IntPermutationsGenerator lowIndicesPermutator, upperIndicesPermutator;
//        Sum sum = new Sum();
//        if (upperCount != 0 && lowerCount != 0) {
//            lowIndicesPermutator = new IntPermutationsGenerator(lowerCount);
//            while (lowIndicesPermutator.hasNext()) {
//                int[] lowerPermutation = lowIndicesPermutator.next().clone();
//                for (int i = 0; i < lowerCount; ++i)
//                    lowerPermutation[i] = lowerPermutation[i] + upperCount;
//                upperIndicesPermutator = new IntPermutationsGenerator(upperCount);
//                UPPER:
//                while (upperIndicesPermutator.hasNext()) {
//                    int[] upperPermutation = upperIndicesPermutator.next();
//                    sum.add(permute(tensor, upperPermutation, lowerPermutation));
//                }
//            }
//        } else if (upperCount == 0) {
//            lowIndicesPermutator = new IntPermutationsGenerator(lowerCount);
//            while (lowIndicesPermutator.hasNext()) {
//                int[] lowerPermutation = lowIndicesPermutator.next();
//                sum.add(permute(tensor, new int[0], lowerPermutation));
//            }
//        } else if (lowerCount == 0) {
//            upperIndicesPermutator = new IntPermutationsGenerator(upperCount);
//            while (upperIndicesPermutator.hasNext()) {
//                int[] upperPermutation = upperIndicesPermutator.next();
//                sum.add(permute(tensor, upperPermutation, new int[0]));
//            }
//        }
//        TensorNumber factor = factor(lowerCount);
//        factor.multiply(factor(upperCount));
//        return new Product(factor, sum);
//    }
//
//    private static TensorNumber factor(int num) {
//        int factor = 1;
//        for (int i = 1; i <= num; ++i)
//            factor *= i;
//        return new TensorNumber(new ComplexElement(factor, 0));
//    }
//
//    private Tensor permute(Tensor tensor, int[] lowerPermutation, int[] upperPermutation) {
//        //creating resulting permutation upper indices are first,
//        //because initial indices are sorted
//        int[] permutation = new int[lowerCount + upperCount];
//        System.arraycopy(upperPermutation, 0, permutation, 0, upperCount);
//        System.arraycopy(lowerPermutation, 0, permutation, upperCount, lowerCount);
//
//
//
//        //processing new indices from permutation
//        final int[] newIndicesNames = new int[indicesNames.length];
//        for (int i = 0; i < indicesNames.length; ++i)
//            newIndicesNames[i] = indicesNames[permutation[i]];
//
//        //processing new tensor
//        IndexMappingImpl im = new IndexMappingImpl(new int[0], indicesNames, newIndicesNames);
//        Tensor current = tensor.clone();
//        ApplyIndexMappingTransformation.INSTANCE.perform(current, im);
//        return current;
//    }
}
