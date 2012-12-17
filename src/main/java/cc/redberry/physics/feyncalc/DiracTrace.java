package cc.redberry.physics.feyncalc;

import cc.redberry.core.indexgenerator.IndexGenerator;
import cc.redberry.core.indices.*;
import cc.redberry.core.number.Complex;
import cc.redberry.core.tensor.*;
import cc.redberry.core.tensor.iterator.TensorLastIterator;
import cc.redberry.core.transformations.ContractIndices;
import cc.redberry.core.transformations.Transformation;
import cc.redberry.core.transformations.expand.Expand;
import cc.redberry.core.transformations.expand.ExpandAll;
import cc.redberry.core.transformations.substitutions.Substitution;

import java.rmi.server.UID;
import java.util.HashMap;

import static cc.redberry.core.indices.IndicesFactory.createSimple;
import static cc.redberry.core.indices.IndicesUtils.setType;
import static cc.redberry.core.tensor.Tensors.*;
import static cc.redberry.core.tensor.Tensors.simpleTensor;

/**
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
public class DiracTrace implements Transformation {
    private final int gammaName, gamma5Name, leviCivitaName;
    private final IndexType metricType, matrixType;
    private final LeviCivitaSimplify LeviCivitaSimplify;

    public DiracTrace() {
        this(parseSimple("G^a'_{a b'}"), parseSimple("G5^a'_b'"), parseSimple("e_{abcd}"));
    }

    public DiracTrace(SimpleTensor gammaMatrix, SimpleTensor gamma5, SimpleTensor leviCivita) {
        //todo check correct input
        this.gammaName = gammaMatrix.getName();
        this.gamma5Name = gamma5.getName();
        this.leviCivitaName = leviCivita.getName();
        IndexType[] types = TraceUtils.extractTypesFromMatrix(gammaMatrix);
        this.metricType = types[0];
        this.matrixType = types[1];
        this.LeviCivitaSimplify = new LeviCivitaSimplify(leviCivita);
    }

    @Override
    public Tensor transform(Tensor tensor) {
        //todo check for contains gammas
        tensor = ExpandAll.expandAll(tensor, ContractIndices.ContractIndices);
        tensor = ContractIndices.contract(tensor);
        TensorLastIterator iterator = new TensorLastIterator(tensor);
        Tensor current;
        while ((current = iterator.next()) != null) {
            if (isGammaOrGamma5(current)
                    && current.getIndices().getFree().size(matrixType) == 0) {
                iterator.set(Complex.ZERO);
            } else if (current instanceof Product) {
                if (current.getIndices().getFree().size(matrixType) != 0)
                    continue;

                int[] gg = calculateGammasInProduct(current);
                int gammasCount = gg[0];
                int gamma5Count = gg[1];

                if (gammasCount == 0) {
                    if (gamma5Count % 2 == 1)
                        iterator.set(Complex.ZERO);
                    else
                        iterator.set(Complex.FOUR);
                    continue;
                }

                if (gamma5Count == 0 && gammasCount % 2 == 1) {
                    iterator.set(Complex.ZERO);
                    continue;
                }

                if (gamma5Count % 2 == 1 && gammasCount % 2 == 1) {
                    iterator.set(Complex.ZERO);
                    continue;
                }

                if (gamma5Count == 0)
                    current = traceWithout5(current, gammasCount);
                else
                    current = trace5((Product) current, gammasCount, gamma5Count);


                //todo d_i^i = 4
                iterator.set(current);
            }
        }

        return iterator.result();
    }

    //first is gamma, second is gamma5
    int[] calculateGammasInProduct(Tensor product) {
        int[] r = new int[2];
        for (Tensor t : product) {
            if (t instanceof SimpleTensor) {
                if (((SimpleTensor) t).getName() == gammaName)
                    ++r[0];
                if (((SimpleTensor) t).getName() == gamma5Name)
                    ++r[1];
            }
        }
        return r;
    }

    /*
     * Trace without gamma5
     */

    private Tensor traceWithout5(Tensor tensor, int gammasCount) {
        Tensor t = createGammaSubstitution(gammasCount).transform(tensor);
        t = Expand.expand(t, ContractIndices.ContractIndices);
        t = ContractIndices.contract(t);
        t = parseExpression("d_a^a=4").transform(t);
        return t;
    }

    //cached trace identities
    private static final HashMap<Key, Substitution> cache = new HashMap<>();

    private static class Key {
        final int gammaName, gammasCount;
        final IndexType metricType, matrixType;

        private Key(int gammaName, int gammasCount, IndexType metricType, IndexType matrixType) {
            this.gammaName = gammaName;
            this.gammasCount = gammasCount;
            this.metricType = metricType;
            this.matrixType = matrixType;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o)
                return true;
            if (o == null || o.getClass() != Key.class)
                return false;
            Key k = (Key) o;
            return gammasCount == k.gammasCount &&
                    gammaName == k.gammaName
                    && metricType == k.metricType
                    && matrixType == k.matrixType;
        }

        @Override
        public int hashCode() {
            int result = gammaName;
            result = 31 * result + metricType.hashCode();
            result = 31 * result + matrixType.hashCode();
            return result;
        }
    }

    private Substitution createGammaSubstitution(int gammasCount) {
        Key key = new Key(gammaName, gammasCount, metricType, matrixType);
        Substitution sub = cache.get(key);
        if (sub != null)
            return sub;
        Tensor[] gammas = new Tensor[gammasCount];
        IndexGenerator generator = new IndexGenerator();
        int firstUpper, u = firstUpper = generator.generate(matrixType), i;
        for (i = 0; i < gammasCount; ++i) {
            gammas[i] = Tensors.simpleTensor(gammaName,
                    createSimple(null,
                            u | 0x80000000,
                            i == gammasCount - 1 ? firstUpper : (u = generator.generate(matrixType)),
                            generator.generate(metricType)));

        }
        sub = new Substitution(Tensors.multiply(gammas), traceOfArray(gammas, metricType));
        cache.put(key, sub);
        return sub;
    }

    static Tensor traceOfArray(Tensor[] product, IndexType metricType) {
        if (product.length == 1)
            return Complex.ZERO;
        if (product.length == 2)
            return multiply(Complex.FOUR,
                    createMetricOrKronecker(product[0].getIndices().get(metricType, 0),
                            product[1].getIndices().get(metricType, 0)));
        if (product.length % 2 != 0)
            return Complex.ZERO;
        SumBuilder sb = new SumBuilder();
        Tensor temp;
        for (int i = 0; i < product.length - 1; ++i) {
            temp = multiply(Complex.TWO,
                    createMetricOrKronecker(product[i].getIndices().get(metricType, 0),
                            product[i + 1].getIndices().get(metricType, 0)),
                    traceOfArray(subArray(product, i, i + 1), metricType));
            if (i % 2 != 0)
                temp = negate(temp);
            sb.put(temp);
            swap(product, i, i + 1);
        }
        return multiply(Complex.ONE_HALF, sb.build());
    }

    /*
     * Trace with gamma5
     */

    private Tensor trace5(Product product, int gammasCount, int gamma5Count) {
        if (gamma5Count == 0)
            return traceWithout5(product, gammasCount);
        if (gamma5Count % 2 == 1 && (gammasCount < 4 || gammasCount % 2 == 1))
            return Complex.ZERO;

        //extracting gamma matrices
        ProductBuilder gammasBuilder = new ProductBuilder(0, gammasCount + gamma5Count);
        int i;
        Tensor t, nonGammaPart = product;
        for (i = product.size() - 1; i >= 0; --i) {
            t = product.get(i);
            if (!isGammaOrGamma5(t))
                continue;

            //removing gamma from product
            if (nonGammaPart instanceof Product)
                nonGammaPart = ((Product) nonGammaPart).remove(i);
            else nonGammaPart = Complex.ONE;

            gammasBuilder.put(t);
        }

        Tensor gammasPart = gammasBuilder.build();
        Substitution[] gamma5 = gamma5Substitution();
        //eliminating multiple gamma5s
        gammasPart = eliminateGamma5((Product) gammasPart, gamma5[0]);

        //even number of gamma5s
        //then all of them were eliminated
        if (gamma5Count % 2 == 0) {
            t = traceWithout5(gammasPart, gammasCount);
            t = multiply(nonGammaPart, t);
            return t;
        }

        //main routine
        t = gamma5[1].transform(gammasPart);
        t = Expand.expand(t);
        t = ContractIndices.contract(t);
        if (t instanceof Sum) {
            SumBuilder sb = new SumBuilder();
            for (Tensor tt : t) {
                int[] gg = calculateGammasInProduct(tt);
                sb.put(trace5((Product) tt, gg[0], gg[1]));
            }
            t = sb.build();
        } else if (t instanceof Product) {
            int[] gg = calculateGammasInProduct(t);
            if (gg[0] != 0 || gg[1] != 0)
                t = trace5((Product) t, gg[0], gg[1]);
        }
        t = multiply(nonGammaPart, t);
        t = Expand.expand(t, ContractIndices.ContractIndices);
        t = ContractIndices.contract(t);
        t = LeviCivitaSimplify.transform(t);
        t = parseExpression("d_a^a=4").transform(t);
        return t;
    }

    private Tensor eliminateGamma5(Product tensor, Substitution sub) {
        //find and remove temporary any non gamma5 matrix:
        Tensor pivot = null;
        Tensor t;
        for (int i = tensor.size() - 1; i >= 0; --i) {
            t = tensor.get(i);
            if (t instanceof SimpleTensor && ((SimpleTensor) t).getName() == gammaName) {
                pivot = t;
                tensor = (Product) tensor.remove(i);
                break;
            }
        }
        Tensor temp = tensor;
        t = tensor;
        do {
            temp = t;
            t = sub.transform(temp);
            t = ContractIndices.contract(t);
        } while (temp != t);

        return multiply(t, pivot);
    }

    private static HashMap<Key5, Substitution[]> cache5 = new HashMap<>();

    private static final class Key5 extends Key {
        final int gamma5Name, leviCivita;

        private Key5(int gammaName, int gammasCount, IndexType metricType, IndexType matrixType, int gamma5Name, int leviCivita) {
            super(gammaName, gammasCount, metricType, matrixType);
            this.gamma5Name = gamma5Name;
            this.leviCivita = leviCivita;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o)
                return true;
            if (o == null || o.getClass() != Key5.class)
                return false;
            Key5 k = (Key5) o;
            return super.equals(o)
                    && gamma5Name == k.gamma5Name
                    && leviCivita == k.leviCivita;
        }

        @Override
        public int hashCode() {
            int result = super.hashCode();
            result = 31 * result + gamma5Name;
            return result;
        }
    }

    //first is swapping, second is trace of 4x1 and Chiholm-Kahane identitie
    private Substitution[] gamma5Substitution() {
        Key5 key = new Key5(gammaName, 3, metricType, matrixType, gamma5Name, leviCivitaName);
        Substitution[] subs = cache5.get(key);
        if (subs != null)
            return subs;
        subs = new Substitution[2];
        Expression[] expressions = new Expression[3];

        //G5*G_m = -G_m*G5
        expressions[0] = expression(
                multiply(
                        simpleTensor(gamma5Name, createSimple(null,
                                setType(matrixType, 0x80000000),
                                setType(matrixType, 1))),
                        simpleTensor(gammaName, createSimple(null,
                                setType(metricType, 0),
                                setType(matrixType, 1 | 0x80000000),
                                setType(matrixType, 2))))
                ,
                multiply(
                        Complex.MINUS_ONE,
                        simpleTensor(gammaName, createSimple(null,
                                setType(metricType, 0),
                                setType(matrixType, 0x80000000),
                                setType(matrixType, 1))),
                        simpleTensor(gamma5Name, createSimple(null,
                                setType(matrixType, 1 | 0x80000000),
                                setType(matrixType, 2))))
        );

        //G5*G5 = 1
        expressions[1] = expression(
                multiply(
                        simpleTensor(gamma5Name, createSimple(null,
                                setType(matrixType, 0x80000000),
                                setType(matrixType, 1))),
                        simpleTensor(gamma5Name, createSimple(null,
                                setType(matrixType, 1 | 0x80000000),
                                setType(matrixType, 2))))
                , createKronecker(setType(matrixType, 0x80000000), setType(matrixType, 2)));

        //Tr[G5] = 0
        expressions[2] = expression(simpleTensor(gamma5Name, createSimple(null,
                setType(matrixType, 0x80000000 | 0),
                setType(matrixType, 0))),
                Complex.ZERO);

        subs[0] = new Substitution(expressions);

        expressions = new Expression[2];

        //Tr[G_a*G_b*G_c*G_d*G5] = -4*I*e_abcd
        ProductBuilder pb = new ProductBuilder(0, 5);
        for (int i = 0; i < 4; ++i) {
            pb.put(simpleTensor(gammaName,
                    createSimple(null,
                            setType(metricType, i),
                            setType(matrixType, 0x80000000 | i),
                            setType(matrixType, i + 1))));
        }
        pb.put(simpleTensor(gamma5Name,
                createSimple(null,
                        setType(matrixType, 0x80000000 | 4),
                        setType(matrixType, 0))));
        expressions[0] = expression(pb.build(),
                multiply(Complex.MINUS_ONE, Complex.FOUR, Complex.IMAGE_ONE,
                        simpleTensor(leviCivitaName, createSimple(null, 0, 1, 2, 3))));

        //Chiholm-Kahane identitie:
        //G_a*G_b*G_c = g_ab*G_c-g_ac*G_b+g_bc*G_a-I*e_abcd*G5*G^d
        Tensor lhs = multiply(
                simpleTensor(gammaName, createSimple(null,
                        setType(metricType, 0),
                        setType(matrixType, 0x80000000),
                        setType(matrixType, 1))),
                simpleTensor(gammaName, createSimple(null,
                        setType(metricType, 1),
                        setType(matrixType, 1 | 0x80000000),
                        setType(matrixType, 2))),
                simpleTensor(gammaName, createSimple(null,
                        setType(metricType, 2),
                        setType(matrixType, 2 | 0x80000000),
                        setType(matrixType, 3)))
        );


        SumBuilder rhs = new SumBuilder();
        Tensor temp;
        for (int i = 0; i < 3; ++i) {
            temp = multiply(simpleTensor(gammaName,
                    createSimple(null,
                            setType(metricType, i),
                            setType(matrixType, 0x80000000),
                            setType(matrixType, 3))),
                    i == 0 ? createMetric(setType(metricType, 1), setType(metricType, 2))
                            : i == 1 ? createMetric(setType(metricType, 0), setType(metricType, 2))
                            : createMetric(setType(metricType, 0), setType(metricType, 1)));
            rhs.put(i == 1 ? negate(temp) : temp);
        }

        temp = multiply(Complex.MINUS_ONE, Complex.IMAGE_ONE,
                simpleTensor(leviCivitaName, createSimple(null,
                        setType(metricType, 0),
                        setType(metricType, 1),
                        setType(metricType, 2),
                        setType(metricType, 3))),
                simpleTensor(gamma5Name, createSimple(null,
                        setType(matrixType, 0x80000000),
                        setType(matrixType, 1))),
                simpleTensor(gammaName, createSimple(null,
                        setType(metricType, 0x80000000 | 3),
                        setType(matrixType, 1 | 0x80000000),
                        setType(matrixType, 3))));
        rhs.put(temp);
        expressions[1] = expression(lhs, rhs.build());
        subs[1] = new Substitution(expressions);
        cache5.put(key, subs);
        return subs;
    }

    private boolean isGammaOrGamma5(Tensor tensor) {
        return tensor instanceof SimpleTensor
                && (((SimpleTensor) tensor).getName() == gammaName || ((SimpleTensor) tensor).getName() == gamma5Name);
    }

    static Tensor traceOfArray5(Tensor[] product, IndexType metricType, int epsilonName) {
        if (product.length < 4 || product.length % 2 == 1)
            return Complex.ZERO;
        if (product.length == 4) {
            int[] indices = new int[4];
            for (int i = 0; i < 4; ++i)
                indices[i] = product[i].getIndices().get(metricType, 0);
            return multiply(Complex.FOUR.negate(), Complex.IMAGE_ONE, simpleTensor(epsilonName, createSimple(null, indices)));
        }
        return null;
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

    public static Tensor trace(Tensor tensor) {
        return new DiracTrace().transform(tensor);
    }
}
