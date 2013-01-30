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
import cc.redberry.core.graph.GraphType;
import cc.redberry.core.graph.PrimitiveSubgraph;
import cc.redberry.core.graph.PrimitiveSubgraphPartition;
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
import cc.redberry.core.transformations.TransformationCollection;
import cc.redberry.core.transformations.expand.ExpandTransformation;
import cc.redberry.core.utils.IntArrayList;
import cc.redberry.core.utils.TensorUtils;

import java.util.ArrayList;

import static cc.redberry.core.tensor.Tensors.multiply;
import static cc.redberry.core.tensor.Tensors.parseExpression;

/**
 * Calculates trace of unitary matrices.
 *
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
    private final Transformation simplifications;

    /**
     * Creates transformation with given definitions.
     *
     * @param unitaryMatrix     unitary matrix
     * @param structureConstant structure constants of SU(N)
     * @param symmetricConstant symmetric constants of SU(N)
     * @param dimension         dimension
     */
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
                        if (!(dimension instanceof Complex))
                            return dimension.toString(OutputFormat.Redberry);
                    default:
                        return old.getName();
                }
            }
        });

        this.unitaryCommutator = (Expression) tokenTransformer.transform(commutatorToken).toTensor();

        //simplifications with SU(N) combinations
        ArrayList<Transformation> unitarySimplifications = new ArrayList<>();
        unitarySimplifications.add(EliminateMetricsTransformation.ELIMINATE_METRICS);
        for (ParseToken substitution : unitarySimplificationsTokens)
            unitarySimplifications.add((Expression) tokenTransformer.transform(substitution).toTensor());
        if (dimension instanceof Complex)
            unitarySimplifications.add(parseExpression("N = " + dimension));

        //all simplifications
        ArrayList<Transformation> simplifications = new ArrayList<>(10);
        simplifications.add(new ExpandTransformation(new TransformationCollection(unitarySimplifications)));
        for (Transformation tr : unitarySimplifications)
            simplifications.add(tr);

        this.simplifications = new TransformationCollection(simplifications);
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
                //extracting trace combinations from product
                Product product = (Product) c;
                int sizeOfIndexless = product.sizeOfIndexlessPart();
                ProductContent productContent = product.getContent();
                PrimitiveSubgraph[] subgraphs
                        = PrimitiveSubgraphPartition.calculatePartition(productContent, matrixType);

                //no unitary matrices in product
                if (subgraphs.length == 0)
                    continue;

                //positions of unitary matrices
                IntArrayList positionsOfMatrices = new IntArrayList();
                //calculated traces
                ProductBuilder calculatedTraces = new ProductBuilder();

                out:
                for (PrimitiveSubgraph subgraph : subgraphs) {
                    //not a trace
                    if (subgraph.getGraphType() != GraphType.Cycle)
                        continue;

                    //positions of unitary matrices
                    int[] partition = subgraph.getPartition();
                    for (int i = partition.length - 1; i >= 0; --i) {
                        partition[i] = sizeOfIndexless + partition[i];
                        //contains not only unitary matrices
                        if (!isUnitaryMatrix(product.get(partition[i]), unitaryMatrix))
                            continue out;
                    }

                    //calculate trace
                    calculatedTraces.put(traceOfProduct(product.select(partition)));
                    positionsOfMatrices.addAll(partition);
                }
                //compiling the result
                c = product.remove(positionsOfMatrices.toArray());
                c = multiply(c, calculatedTraces.build());
                iterator.set(simplifications.transform(c));
            }
        }
        return iterator.result();
    }

    private Tensor traceOfProduct(Tensor tensor) {
        Tensor oldTensor = tensor, newTensor;
        while (true) {
            newTensor = unitaryCommutator.transform(oldTensor);
            newTensor = simplifications.transform(newTensor);
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

        if (dimension instanceof Complex && !TensorUtils.isNaturalNumber(dimension))
            throw new IllegalArgumentException("Non natural dimension.");

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
                    throw new IllegalArgumentException("Different indices metric types: " + unitaryMatrix + " and " + st);
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
        numberOfGeneratorsToken = parser.parse("d^a_a = N**2-1");
        dimensionToken = parser.parse("d^a'_a' = N");

        unitarySimplificationsTokens = new ParseToken[]{
                singleTraceToken, symmetricCombinationToken, aSymmetricCombinationToken,
                symmetrySimplificationToken, numberOfGeneratorsToken, dimensionToken};
    }


}
