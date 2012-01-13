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

import java.util.List;
import redberry.core.number.ComplexElement;
import redberry.core.number.NumberFraction;
import redberry.core.number.RationalElement;
import redberry.core.tensor.Product;
import redberry.core.tensor.Sum;
import redberry.core.tensor.Tensor;
import redberry.core.tensor.TensorNumber;
import redberry.core.tensor.generator.IndexMappingPermutationsGenerator;
import redberry.core.transformation.Transformation;

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
//    private final Indexes indexes;
//    private final int lowerCount;
//    private final int upperCount;
//    private final int[] indexesNames;
//
//    public Symmetrize(Indexes indexes) {
//        this.indexes = indexes.getFreeIndexes();
//        lowerCount = indexes.getLower().size();
//        upperCount = indexes.getUpper().size();
//        indexesNames = new int[this.indexes.size()];
//        for (int i = 0; i < indexesNames.length; ++i)
//            indexesNames[i] = IndexesUtils.getNameWithType(this.indexes.get(i));
//    }
//
//    public Symmetrize(String indexes) {
//        this(ParserIndexes.parse(indexes));
//    }
//
//    @Override
//    public Tensor transform(Tensor tensor) {
//        if (!tensor.getIndexes().getFreeIndexes().equalsIgnoreOrder(indexes))
//            throw new IllegalArgumentException();
//
//        IntPermutationsGenerator lowIndexesPermutator, upperIndexesPermutator;
//        Sum sum = new Sum();
//        if (upperCount != 0 && lowerCount != 0) {
//            lowIndexesPermutator = new IntPermutationsGenerator(lowerCount);
//            while (lowIndexesPermutator.hasNext()) {
//                int[] lowerPermutation = lowIndexesPermutator.next().clone();
//                for (int i = 0; i < lowerCount; ++i)
//                    lowerPermutation[i] = lowerPermutation[i] + upperCount;
//                upperIndexesPermutator = new IntPermutationsGenerator(upperCount);
//                UPPER:
//                while (upperIndexesPermutator.hasNext()) {
//                    int[] upperPermutation = upperIndexesPermutator.next();
//                    sum.add(permute(tensor, upperPermutation, lowerPermutation));
//                }
//            }
//        } else if (upperCount == 0) {
//            lowIndexesPermutator = new IntPermutationsGenerator(lowerCount);
//            while (lowIndexesPermutator.hasNext()) {
//                int[] lowerPermutation = lowIndexesPermutator.next();
//                sum.add(permute(tensor, new int[0], lowerPermutation));
//            }
//        } else if (lowerCount == 0) {
//            upperIndexesPermutator = new IntPermutationsGenerator(upperCount);
//            while (upperIndexesPermutator.hasNext()) {
//                int[] upperPermutation = upperIndexesPermutator.next();
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
//        //creating resulting permutation upper indexes are first,
//        //because initial indexes are sorted
//        int[] permutation = new int[lowerCount + upperCount];
//        System.arraycopy(upperPermutation, 0, permutation, 0, upperCount);
//        System.arraycopy(lowerPermutation, 0, permutation, upperCount, lowerCount);
//
//
//
//        //processing new indexes from permutation
//        final int[] newIndexesNames = new int[indexesNames.length];
//        for (int i = 0; i < indexesNames.length; ++i)
//            newIndexesNames[i] = indexesNames[permutation[i]];
//
//        //processing new tensor
//        IndexMappingImpl im = new IndexMappingImpl(new int[0], indexesNames, newIndexesNames);
//        Tensor current = tensor.clone();
//        ApplyIndexMappingTransformation.INSTANCE.perform(current, im);
//        return current;
//    }
}
