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

import cc.redberry.core.context.CC;
import cc.redberry.core.context.NameAndStructureOfIndices;
import cc.redberry.core.context.OutputFormat;
import cc.redberry.core.indices.IndexType;
import cc.redberry.core.indices.IndicesUtils;
import cc.redberry.core.number.Complex;
import cc.redberry.core.parser.ParseToken;
import cc.redberry.core.parser.Parser;
import cc.redberry.core.parser.preprocessor.ChangeIndicesTypesAndTensorNames;
import cc.redberry.core.parser.preprocessor.TypesAndNamesTransformer;
import cc.redberry.core.tensor.*;
import cc.redberry.core.tensor.iterator.FromChildToParentIterator;
import cc.redberry.core.transformations.EliminateMetricsTransformation;
import cc.redberry.core.transformations.Transformation;
import cc.redberry.core.transformations.expand.ExpandTransformation;
import cc.redberry.core.utils.TensorUtils;

import static cc.redberry.core.tensor.Tensors.multiply;

/**
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
public final class UnitaryTrace implements Transformation {
    /*
     * Defaults
     */
    private static final String unitaryMatrixName = "T";
    private static final String structureConstantName = "F";
    private static final String symmetricConstantName = "D";
    private static final String dimensionName = "N";

    private final int unitaryMatrix;
    private final IndexType matrixType;

    private final Expression unitaryCommutator;
    private final Transformation[] simplifications;

    public UnitaryTrace(final SimpleTensor unitaryMatrix,
                        final SimpleTensor structureConstant,
                        final SimpleTensor symmetricConstant,
                        final Tensor dimension) {
        check(unitaryMatrix, structureConstant, symmetricConstant, dimension);
        this.unitaryMatrix = unitaryMatrix.getName();
        final IndexType[] types = TraceUtils.extractTypesFromMatrix(unitaryMatrix);
        this.matrixType = types[1];

        ChangeIndicesTypesAndTensorNames tokenTransformer = new ChangeIndicesTypesAndTensorNames(new TypesAndNamesTransformer() {
            @Override
            public IndexType newType(IndexType oldType, NameAndStructureOfIndices old) {
                if (oldType == IndexType.LatinLower)
                    return types[0];
                if (oldType == IndexType.Matrix1)
                    return types[1];
                return oldType;
            }

            @Override
            public String newName(NameAndStructureOfIndices old) {
                switch (old.getName()) {
                    case unitaryMatrixName:
                        return unitaryMatrix.getStringName();
                    case structureConstantName:
                        return structureConstant.getStringName();
                    case symmetricConstantName:
                        return symmetricConstant.getStringName();
                    case dimensionName:
                        return dimension.toString(OutputFormat.Redberry);
                    default:
                        return old.getName();
                }
            }
        });

        this.unitaryCommutator = (Expression) tokenTransformer.transform(commutatorToken).toTensor();
        this.simplifications = new Transformation[1 + unitarySimplificationsTokens.length];

        int i = 0;
        this.simplifications[i++] = EliminateMetricsTransformation.ELIMINATE_METRICS;
        for (ParseToken substitution : unitarySimplificationsTokens)
            this.simplifications[i++] = (Expression) tokenTransformer.transform(substitution).toTensor();
    }

    @Override
    public Tensor transform(Tensor t) {
        FromChildToParentIterator iterator = new FromChildToParentIterator(t);
        Tensor c;
        while ((c = iterator.next()) != null) {
            if (c instanceof SimpleTensor) {
                // Tr[T_a] = 0
                if (((SimpleTensor) c).getName() == unitaryMatrix && c.getIndices().getOfType(matrixType).getFree().size() == 0)
                    iterator.set(Complex.ZERO);
            } else if (c instanceof Product) {
                //selecting unitary matrices from product

                //subproduct of matrices
                TensorBuilder matrices = new ProductBuilder();
                boolean noUnitaryMatrices = true;
                //remainder (non matrix part)
                Tensor remainder = c;
                for (int i = c.size() - 1; i >= 0; --i) {
                    if (isUnitaryMatrix(c.get(i), unitaryMatrix)) {
                        noUnitaryMatrices = false;
                        matrices.put(c.get(i));
                        if (remainder instanceof Product)
                            remainder = ((Product) remainder).remove(i);
                        else {
                            assert i == 0;
                            remainder = Complex.ONE;
                        }
                    }
                }
                if (noUnitaryMatrices)
                    continue;

                //compiling the result
                Tensor temp = multiply(remainder, traceOfProduct(matrices.build()));
                temp = ExpandTransformation.expand(temp, EliminateMetricsTransformation.ELIMINATE_METRICS);
                temp = EliminateMetricsTransformation.eliminate(temp);
                iterator.set(temp);
            }
        }
        return iterator.result();
    }

    private Tensor traceOfProduct(Tensor tensor) {
        Tensor oldTensor = tensor, newTensor;
        while (true) {
            newTensor = unitaryCommutator.transform(oldTensor);
            newTensor = ExpandTransformation.expand(newTensor, simplifications);
            for (Transformation tr : simplifications)
                newTensor = tr.transform(newTensor);
            if (newTensor == oldTensor)
                break;
            oldTensor = newTensor;
        }
        return newTensor;
    }

    private static final boolean isUnitaryMatrix(Tensor tensor, int unitaryMatrix) {
        return tensor instanceof SimpleTensor && ((SimpleTensor) tensor).getName() == unitaryMatrix;
    }

    private static void check(final SimpleTensor unitaryMatrix,
                              final SimpleTensor structureConstant,
                              final SimpleTensor symmetricConstant,
                              final Tensor dimension) {
        if (unitaryMatrix.getIndices().size() != 3)
            throw new IllegalArgumentException("Not a unitary matrix: " + unitaryMatrix);
        IndexType[] types = TraceUtils.extractTypesFromMatrix(unitaryMatrix);
        IndexType metricType = types[0];
        if (!TensorUtils.isScalar(dimension))
            throw new IllegalArgumentException("Non scalar dimension.");
        if (structureConstant.getName() == symmetricConstant.getName())
            throw new IllegalArgumentException("Structure and symmetric constants have same names.");
        SimpleTensor[] ss = {structureConstant, symmetricConstant};
        for (SimpleTensor st : ss) {
            if (st.getIndices().size() != 3)
                throw new IllegalArgumentException("Illegal input for SU(N) constants: " + st);
            for (int i = 0; i < 3; ++i)
                if (IndicesUtils.getTypeEnum(st.getIndices().get(i)) != metricType)
                    throw new IllegalArgumentException("Different indices types: " + unitaryMatrix + " and " + st);
        }
    }

    /*
     * Substitutions
     */

    private static final Parser parser;
    /**
     * T_a*T_b  = 1/2N g_ab + I/2*f_abc*T^c + 1/2*d_abc*T^c
     */
    private static final ParseToken commutatorToken;
    /**
     * Tr[T_a] = 0
     */
    private static final ParseToken singleTraceToken;

    /**
     * d_apq*d_b^pq = (N**2 - 4)/N * g_ab
     */
    private static final ParseToken symmetricCombinationToken;

    /**
     * f_apq*f_b^pq = N * g_ab
     */
    private static final ParseToken aSymmetricCombinationToken;

    /**
     * f_apq*d_b^pq = 0
     */
    private static final ParseToken symmetrySimplificationToken;

    /**
     * d^a_a = N*(N-1)/2
     */
    private static final ParseToken numberOfGeneratorsToken;

    /**
     * Tr[1] = N
     */
    private static final ParseToken dimensionToken;

    private static final ParseToken[] unitarySimplificationsTokens;

    static {
        parser = CC.current().getParseManager().getParser();

        commutatorToken = parser.parse("T_a^a'_c'*T_b^c'_b' = 1/(2*N)*g_ab*d^a'_b' + I/2*F_abc*T^ca'_b' + 1/2*D_abc*T^ca'_b'");
        singleTraceToken = parser.parse("T_a^a'_a' = 0");
        symmetricCombinationToken = parser.parse("D_apq*D_b^pq = (N**2 - 4)/N * g_ab");
        aSymmetricCombinationToken = parser.parse("F_apq*F_b^pq = N * g_ab");
        symmetrySimplificationToken = parser.parse("F_apq*D_b^pq = 0");
        numberOfGeneratorsToken = parser.parse("d^a_a = N*(N-1)/2");
        dimensionToken = parser.parse("d^a'_a' = N");

        unitarySimplificationsTokens = new ParseToken[]{
                singleTraceToken, symmetricCombinationToken, aSymmetricCombinationToken,
                symmetrySimplificationToken, numberOfGeneratorsToken, dimensionToken};
    }


}
