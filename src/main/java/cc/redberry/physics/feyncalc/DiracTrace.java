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
import cc.redberry.core.graph.GraphType;
import cc.redberry.core.graph.PrimitiveSubgraph;
import cc.redberry.core.graph.PrimitiveSubgraphPartition;
import cc.redberry.core.indices.IndexType;
import cc.redberry.core.indices.IndicesFactory;
import cc.redberry.core.number.Complex;
import cc.redberry.core.parser.ParseToken;
import cc.redberry.core.parser.ParseTokenTransformer;
import cc.redberry.core.parser.ParseUtils;
import cc.redberry.core.parser.Parser;
import cc.redberry.core.parser.preprocessor.ChangeIndicesTypesAndTensorNames;
import cc.redberry.core.parser.preprocessor.TypesAndNamesTransformer;
import cc.redberry.core.tensor.*;
import cc.redberry.core.tensor.iterator.FromChildToParentIterator;
import cc.redberry.core.transformations.ExpandAndEliminateTransformation;
import cc.redberry.core.transformations.Transformation;
import cc.redberry.core.transformations.TransformationCollection;
import cc.redberry.core.utils.IntArrayList;
import gnu.trove.map.hash.TIntObjectHashMap;

import java.util.ArrayList;

import static cc.redberry.core.indices.IndicesFactory.createSimple;
import static cc.redberry.core.indices.IndicesUtils.*;
import static cc.redberry.core.tensor.Tensors.*;

/**
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
public class DiracTrace implements Transformation {
    /*
     * Defaults
     */
    private static final String gammaMatrixStringName = "G";
    private static final String gamma5StringName = "G5";
    private static final String leviCivitaStringName = "eps";

    private final int gammaName, gamma5Name, leviCivitaName;
    private final IndexType metricType, matrixType;
    private final LeviCivitaSimplify LeviCivitaSimplify;
    private final boolean minkovskiSpace;

    private final ParseTokenTransformer tokenTransformer;

    public DiracTrace(final SimpleTensor gammaMatrix,
                      final SimpleTensor gamma5,
                      final SimpleTensor leviCivita,
                      final boolean minkovskiSpace) {
        //todo check correct input
        this.minkovskiSpace = minkovskiSpace;
        this.gammaName = gammaMatrix.getName();
        this.gamma5Name = gamma5.getName();
        this.leviCivitaName = leviCivita.getName();
        final IndexType[] types = TraceUtils.extractTypesFromMatrix(gammaMatrix);
        this.metricType = types[0];
        this.matrixType = types[1];

        tokenTransformer = new ChangeIndicesTypesAndTensorNames(new TypesAndNamesTransformer() {
            @Override
            public IndexType newType(IndexType oldType, NameAndStructureOfIndices oldDescriptor) {
                switch (oldType) {
                    case LatinLower:
                        return metricType;
                    case Matrix1:
                        return matrixType;
                }
                return oldType;
            }

            @Override
            public String newName(NameAndStructureOfIndices oldDescriptor) {
                switch (oldDescriptor.getName()) {
                    case gammaMatrixStringName:
                        return gammaMatrix.getStringName();
                    case gamma5StringName:
                        return gamma5StringName;
                    case leviCivitaStringName:
                        return leviCivita.getStringName();
                    default:
                        return oldDescriptor.getName();
                }
            }
        });

        this.LeviCivitaSimplify = new LeviCivitaSimplify(leviCivita, minkovskiSpace);

    }

    @Override
    public Tensor transform(Tensor tensor) {
        //todo check for contains gammas
        tensor = ExpandAndEliminateTransformation.expandAndEliminate(tensor);
        FromChildToParentIterator iterator = new FromChildToParentIterator(tensor);
        Tensor current;
        out:
        while ((current = iterator.next()) != null) {
            if (isGammaOrGamma5(current)
                    && current.getIndices().getFree().size(matrixType) == 0) {
                iterator.set(Complex.ZERO);
            } else if (current instanceof Product) {
                if (current.getIndices().getFree().size(matrixType) != 0)
                    continue;
                //selecting unitary matrices from product
                //extracting trace combinations from product
                Product product = (Product) current;
                //positions of matrices
                IntArrayList positionsOfMatrices = new IntArrayList();
                int sizeOfIndexless = product.sizeOfIndexlessPart();
                ProductContent pc = product.getContent();
                PrimitiveSubgraph[] partition
                        = PrimitiveSubgraphPartition.calculatePartition(pc, matrixType);

                //calculated traces
                ProductBuilder traces = new ProductBuilder();

                traces:
                for (PrimitiveSubgraph subgraph : partition) {
                    if (subgraph.getGraphType() != GraphType.Cycle)
                        continue;

                    int numberOfGammas = 0, numberOfGamma5s = 0;

                    Tensor gamma;
                    //actual positions in current
                    int[] positions = subgraph.getPartition();
                    assert positions.length > 1;

                    for (int i = positions.length - 1; i >= 0; --i) {
                        positions[i] = positions[i] + sizeOfIndexless;
                        gamma = product.get(positions[i]);
                        if (gamma instanceof SimpleTensor) {
                            if (((SimpleTensor) gamma).getName() == gammaName)
                                ++numberOfGammas;
                            else if (((SimpleTensor) gamma).getName() == gamma5Name)
                                ++numberOfGamma5s;
                            else
                                //not a gamma matrix
                                continue traces;
                        } else {
                            //not a gamma matrix
                            continue traces;
                        }
                    }

                    //early terminations
                    if (numberOfGammas % 2 == 1
                            || (numberOfGammas == 2 && numberOfGamma5s % 2 == 1)) {
                        iterator.set(Complex.ZERO);
                        continue out;
                    }
                    if (numberOfGammas == 0 && numberOfGamma5s % 2 == 1) {
                        iterator.set(Complex.ZERO);
                        continue out;
                    }

                    positionsOfMatrices.addAll(positions);
                    if (numberOfGamma5s == 0)
                        traces.put(traceWithout5(product.select(positions), numberOfGammas));
                    else {
                        //early check
                        if (numberOfGammas == 0) {
                            //numberOfGamma5s % 2 == 0
                            traces.put(Complex.FOUR);
                            continue traces;
                        }


                        //eliminating excess products of gamma5s
                        if (numberOfGamma5s > 1) {
                            //take into account odd number of swaps
                            boolean sign = false;
                            //product of gammas as ordered array (will be filled without excess gamma5s)
                            final SimpleTensor[] orderedProduct = new SimpleTensor[numberOfGammas + (numberOfGamma5s % 2 == 0 ? 0 : 1)];
                            int counter = -1;

                            //index of tensor in product content, which is contracted with current gamma5
                            int positionOfPreviousGamma = -2;

                            SimpleTensor currentGamma;
                            for (int positionOfGamma : positions) {
                                currentGamma = (SimpleTensor) product.get(positionOfGamma);
                                if (currentGamma.getName() == gamma5Name) {
                                    //adding one gamma5 if they are odd number
                                    if (positionOfPreviousGamma == -2 && numberOfGamma5s % 2 == 1) {
                                        orderedProduct[++counter] = currentGamma;
                                        positionOfPreviousGamma = -1;
                                        continue;
                                    }
                                    if (positionOfPreviousGamma == -1)
                                        positionOfPreviousGamma = positionOfGamma;
                                    else {
                                        //odd number of swaps
                                        if ((positionOfGamma - positionOfPreviousGamma) % 2 == 1)
                                            sign ^= sign;
                                        positionOfPreviousGamma = positionOfGamma;
                                    }
                                } else
                                    orderedProduct[++counter] = currentGamma;
                            }

                            //fixing new indices contractions
                            int u = 0, l = 0;
                            for (int i = 0; ; ++i) {
                                if (i == orderedProduct.length - 1) {
                                    orderedProduct[i] = setMatrixIndices(orderedProduct[i], u, 0);
                                    break;
                                }
                                orderedProduct[i] = setMatrixIndices(orderedProduct[i], u, ++l);
                                u = l;
                            }

                            Tensor withoutExcessGamma5s = multiply(orderedProduct);

                            if (numberOfGamma5s % 2 == 0)
                                withoutExcessGamma5s = traceWithout5(withoutExcessGamma5s, numberOfGammas);
                            else
                                withoutExcessGamma5s = traceWith5(withoutExcessGamma5s, numberOfGammas);

                            if (sign)
                                withoutExcessGamma5s = negate(withoutExcessGamma5s);
                            traces.put(withoutExcessGamma5s);
                        } else
                            traces.put(traceWith5(product.select(positions), numberOfGammas));
                    }

                }

                //final simplifications
                traces.put(product.remove(positionsOfMatrices.toArray()));
                iterator.set(ExpandAndEliminateTransformation.expandAndEliminate(traces.build()));
            }
        }
        return iterator.result();
    }

    private SimpleTensor setMatrixIndices(SimpleTensor gamma, int matrixUpper, int matrixLower) {
        int[] indices = gamma.getIndices().getAllIndices().copy();
        for (int i = indices.length - 1; i >= 0; --i)
            if (!CC.isMetric(getType(indices[i]))) {
                indices[i] = getState(indices[i]) ?
                        createIndex(matrixUpper, getType(indices[i]), getState(indices[i]))
                        : setType(getType(indices[i]), matrixLower);
            }
        return simpleTensor(gamma.getName(), IndicesFactory.createSimple(null, indices));
    }


    private boolean isGammaOrGamma5(Tensor tensor) {
        return tensor instanceof SimpleTensor
                && (((SimpleTensor) tensor).getName() == gammaName || ((SimpleTensor) tensor).getName() == gamma5Name);
    }

    /*
     * *********************
     * Trace without gamma5
     * *********************
     */

    /**
     * cached substitutions of traces without 5
     */
    private final TIntObjectHashMap<Expression> cachedSubstitutions = new TIntObjectHashMap<>();

    private Tensor traceWithout5(Tensor tensor, int numberOfGammas) {
        Expression substitution = cachedSubstitutions.get(numberOfGammas);
        if (substitution == null) {
            ParseToken rawSubstitution = createRawGammaSubstitution(numberOfGammas);
            substitution = (Expression) tokenTransformer.transform(rawSubstitution).toTensor();
            cachedSubstitutions.put(numberOfGammas, substitution);
        }

        tensor = substitution.transform(tensor);
        tensor = ExpandAndEliminateTransformation.expandAndEliminate(tensor);
        tensor = parseExpression("d^a_a = 4").transform(tensor);
        return tensor;
    }

    private static ParseToken createRawGammaSubstitution(int numberOfGammas) {
        ParseToken substitution = cachedRawGammaTraces.get(numberOfGammas);
        if (substitution == null) {

            //product of gamma matrices as array
            Tensor[] gammas = new Tensor[numberOfGammas];
            int matrixIndex = setType(IndexType.Matrix1, 0) - 1, metricIndex = -1;
            int firstUpper, u = firstUpper = ++matrixIndex, i;
            for (i = 0; i < numberOfGammas; ++i) {
                gammas[i] = Tensors.simpleTensor(gammaMatrixStringName,
                        createSimple(null,
                                u | 0x80000000,
                                i == numberOfGammas - 1 ? firstUpper : (u = ++matrixIndex),
                                ++metricIndex));

            }
            Expression expression = expression(Tensors.multiply(gammas), traceOfArray(gammas));
            substitution = ParseUtils.tensor2AST(expression);
            cachedRawGammaTraces.put(numberOfGammas, substitution);
        }
        return substitution;
    }

    private static Tensor traceOfArray(Tensor[] product) {
        //calculates trace using recursive algorithm
        if (product.length == 1)
            return Complex.ZERO;
        if (product.length == 2)
            return multiply(Complex.FOUR,
                    createMetricOrKronecker(product[0].getIndices().get(IndexType.LatinLower, 0),
                            product[1].getIndices().get(IndexType.LatinLower, 0)));
        if (product.length % 2 != 0)
            return Complex.ZERO;
        SumBuilder sb = new SumBuilder();
        Tensor temp;
        for (int i = 0; i < product.length - 1; ++i) {
            temp = multiply(Complex.TWO,
                    createMetricOrKronecker(product[i].getIndices().get(IndexType.LatinLower, 0),
                            product[i + 1].getIndices().get(IndexType.LatinLower, 0)),
                    traceOfArray(subArray(product, i, i + 1)));
            if (i % 2 != 0)
                temp = negate(temp);
            sb.put(temp);
            swap(product, i, i + 1);
        }
        return multiply(Complex.ONE_HALF, sb.build());
    }

    private static Tensor[] subArray(Tensor[] array, int a, int b) {
        Tensor[] result = new Tensor[array.length - 2];
        int k = 0;
        for (int i = 0; i < array.length; ++i) {
            if (i == a || i == b)
                continue;
            result[k++] = array[i];
        }
        return result;
    }

    private static void swap(Tensor[] array, int a, int b) {
        Tensor temp = array[a];
        array[a] = array[b];
        array[b] = temp;
    }

    /*
    * Cached parse tokens
    */
    private static final TIntObjectHashMap<ParseToken> cachedRawGammaTraces = new TIntObjectHashMap<>();


    /*
     * *********************
     * Trace with gamma5
     * *********************
     */

    //cached substitutions of traces with 5
    private Transformation cachedSubstitutions5;


    private Tensor traceWith5(Tensor product, int numberOfGammas) {
        if (cachedSubstitutions5 == null) {
            ArrayList<Transformation> substitutions = new ArrayList<>();
            substitutions.add((Expression) tokenTransformer.transform(traceOf4GammasWith5Tokens).toTensor());
            substitutions.add((Expression) tokenTransformer.transform(chiholmKahaneTokens).toTensor());
            cachedSubstitutions5 = new TransformationCollection(substitutions);
        }

        if (numberOfGammas == 4)
            return cachedSubstitutions5.transform(product);

        product = cachedSubstitutions5.transform(product);
        product = ExpandAndEliminateTransformation.expandAndEliminate(product);
        product = cachedSubstitutions5.transform(product);
        product = LeviCivitaSimplify.transform(product);
        return transform(product);
    }


    private static final Parser parser;
    /**
     * Tr[G_a*G_b*G_c*G_d*G5] = -4*I*e_abcd
     */
    private static final ParseToken traceOf4GammasWith5Tokens;
    /**
     * Chiholm-Kahane identitie:
     * G_a*G_b*G_c = g_ab*G_c-g_ac*G_b+g_bc*G_a-I*e_abcd*G5*G^d
     */
    private static final ParseToken chiholmKahaneTokens;

    static {
        parser = CC.current().getParseManager().getParser();
        traceOf4GammasWith5Tokens = parser.parse("G_a^a'_b'*G_b^b'_c'*G_c^c'_d'*G_d^d'_e'*G5^e'_a' = -4*I*eps_abcd");
        chiholmKahaneTokens = parser.parse("G_a^a'_c'*G_b^c'_d'*G_c^d'_b' = g_ab*G_c^a'_b'-g_ac*G_b^a'_b'+g_bc*G_a^a'_b'-I*e_abcd*G5^a'_c'*G^dc'_b'");
    }
}
