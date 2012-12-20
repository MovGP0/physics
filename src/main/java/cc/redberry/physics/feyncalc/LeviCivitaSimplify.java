package cc.redberry.physics.feyncalc;

import cc.redberry.core.combinatorics.Combinatorics;
import cc.redberry.core.combinatorics.Permutation;
import cc.redberry.core.combinatorics.Symmetry;
import cc.redberry.core.combinatorics.symmetries.Symmetries;
import cc.redberry.core.combinatorics.symmetries.SymmetriesFactory;
import cc.redberry.core.indexgenerator.IndexGenerator;
import cc.redberry.core.indexmapping.IndexMappingBuffer;
import cc.redberry.core.indexmapping.IndexMappings;
import cc.redberry.core.indexmapping.MappingsPort;
import cc.redberry.core.indices.IndicesFactory;
import cc.redberry.core.indices.SimpleIndices;
import cc.redberry.core.number.Complex;
import cc.redberry.core.tensor.*;
import cc.redberry.core.tensor.iterator.TensorLastIterator;
import cc.redberry.core.transformations.ContractIndices;
import cc.redberry.core.transformations.Transformation;
import cc.redberry.core.transformations.expand.Expand;
import cc.redberry.core.utils.IntArray;
import cc.redberry.core.utils.IntArrayList;
import cc.redberry.core.utils.TensorUtils;
import gnu.trove.map.hash.TIntObjectHashMap;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import static cc.redberry.core.indices.IndicesUtils.getType;
import static cc.redberry.core.indices.IndicesUtils.inverseIndexState;
import static cc.redberry.core.tensor.FullContractionsStructure.getToTensorIndex;
import static cc.redberry.core.tensor.Tensors.*;

/**
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
public class LeviCivitaSimplify implements Transformation {
    private final SimpleTensor LeviCivita;

    public LeviCivitaSimplify(SimpleTensor leviCivita) {
        checkLeviCivita(leviCivita);
        LeviCivita = leviCivita;
    }

    @Override
    public Tensor transform(Tensor t) {
        return simplifyLeviCivita1(t, LeviCivita);
    }

    public static Tensor simplifyLeviCivita(Tensor tensor, SimpleTensor LeviCivita) {
        checkLeviCivita(LeviCivita);
        return simplifyLeviCivita1(tensor, LeviCivita);
    }

    private static Tensor simplifyLeviCivita1(Tensor tensor, SimpleTensor LeviCivita) {
        TensorLastIterator iterator = new TensorLastIterator(tensor);
        Tensor c;
        while ((c = iterator.next()) != null) {
            if (c instanceof SimpleTensor
                    && ((SimpleTensor) c).getName() == LeviCivita.getName()
                    && c.getIndices().size() != c.getIndices().getFree().size()) {
                iterator.set(Complex.ZERO);
            }
            if (c instanceof Product)
                iterator.set(simplifyProduct(c, LeviCivita));
        }
        return iterator.result();
    }

    private static void checkLeviCivita(SimpleTensor LeviCivita) {
        SimpleIndices indices = LeviCivita.getIndices();
        if (indices.size() <= 1)
            throw new IllegalArgumentException("Levi-Civita cannot be a scalar.");
        byte type = getType(indices.get(0));
        for (int i = 1; i < indices.size(); ++i)
            if (type != getType(indices.get(i)))
                throw new IllegalArgumentException("Levi-Civita have indices with different types.");
    }

    private static Tensor simplifyProduct(Tensor product, SimpleTensor LeviCivita) {
        ProductContent content = ((Product) product).getContent();
        IntArrayList epsPositions = new IntArrayList();
        int i = 0, j = content.size();
        for (i = 0; i < j; ++i) {
            if (isLeviCivita(content.get(i), LeviCivita))
                epsPositions.add(i);
        }
        if (epsPositions.isEmpty())
            return product;
        FullContractionsStructure fs = content.getFullContractionsStructure();

        j = epsPositions.size();
        Set<Tensor> epsComponent = new HashSet<>(LeviCivita.getIndices().size());
        Tensor temp;
        int toIndex, a, b;
        for (i = 0; i < j; ++i) {
            for (long contraction : fs.contractions[epsPositions.get(i)]) {
                toIndex = getToTensorIndex(contraction);
                if (toIndex == -1)
                    continue;
                temp = content.get(toIndex);
                if (isLeviCivita(temp, LeviCivita))
                    continue;
                epsComponent.add(temp);
            }
            if (epsComponent.isEmpty())
                continue;

            temp = multiply(epsComponent.toArray(new Tensor[epsComponent.size()]));
            epsComponent.clear();
            MappingsPort port = IndexMappings.createPort(temp, temp);
            IndexMappingBuffer buffer;
            Symmetry sym;

            IntArrayList nonPermutable = new IntArrayList();
            int[] indices = temp.getIndices().getFree().getAllIndices().copy();
            int[] epsIndices = content.get(epsPositions.get(i)).getIndices().getFree().getAllIndices().copy();

            boolean contract;
            for (b = 0; b < indices.length; ++b) {
                contract = false;
                for (a = 0; a < epsIndices.length; ++a)
                    if (indices[b] == inverseIndexState(epsIndices[a]))
                        contract = true;
                if (!contract)
                    nonPermutable.add(b);
            }
            int[] nonPermutablePositions = nonPermutable.toArray();

            if (indices.length == 1)
                continue;
            Map<IntArray, Boolean> symmetries = getEpsilonSymmetries(indices.length);
            while ((buffer = port.take()) != null) {
                sym = TensorUtils.getSymmetryFromMapping(indices, buffer);
                if (!checkNonPermutingPositions(sym, nonPermutablePositions))
                    continue;

                if (sym.isAntiSymmetry() != symmetries.get(sym.getPermutation()))
                    return Complex.ZERO;
            }

        }
        if (epsPositions.size() == 1)
            return product;

        Expression[] subs = getLeviCivitaSubstitutions(LeviCivita);
        for (Expression exp : subs)
            product = exp.transform(product);
        //todo expand only Levi-Civita sums
        product = ContractIndices.contract(Expand.expand(product, ContractIndices.ContractIndices));
        product = subs[1].transform(product);
        return product;
    }

    private static boolean checkNonPermutingPositions(Permutation permutation, int[] nonPermutablePositions) {
        for (int i : nonPermutablePositions)
            if (permutation.newIndexOf(i) != i)
                return false;
        return true;
    }

    private static boolean isLeviCivita(Tensor tensor, SimpleTensor LeviCivita) {
        return tensor instanceof SimpleTensor && ((SimpleTensor) tensor).getName() == LeviCivita.getName();
    }

    public static Expression[] getLeviCivitaSubstitutions(SimpleTensor eps) {
        Expression[] sub = cachedSubstitutions.get(eps.getName());
        if (sub != null)
            return sub;
        sub = new Expression[2];
        SimpleIndices indices = eps.getIndices();
        int size = indices.size();

        SimpleTensor eps1 = setDiffIndices(eps),
                eps2 = setInversedIndices(setDiffIndices(eps1));
        Tensor lhs = multiply(eps1, eps2);
        SimpleIndices eps1Indices = eps1.getIndices(),
                eps2Indices = eps2.getIndices();
        Tensor[][] matrix = new Tensor[size][size];
        int j;
        for (int i = 0; i < size; ++i)
            for (j = 0; j < size; ++j)
                matrix[i][j] = createKronecker(eps1Indices.get(i), eps2Indices.get(j));

        Tensor rhs = TensorUtils.det(matrix);
        //todo here we assuming, that all indices belongs to pseudo-Euclidean spaces
        if (size % 2 == 0)
            rhs = negate(rhs);
        sub[0] = expression(lhs, rhs);
        int index = eps1Indices.get(0);
        sub[1] = expression(createKronecker(index, inverseIndexState(index)), new Complex(size));
        cachedSubstitutions.put(eps.getName(), sub);
        return sub;
    }

    private static SimpleTensor setDiffIndices(SimpleTensor eps) {
        byte type = getType(eps.getIndices().get(0));
        int count = eps.getIndices().size();
        int[] newIndices = new int[count];
        IndexGenerator generator = new IndexGenerator(eps.getIndices());
        for (int i = 0; i < count; ++i)
            newIndices[i] = generator.generate(type);
        return Tensors.simpleTensor(eps.getName(), IndicesFactory.createSimple(null, newIndices));
    }

    private static SimpleTensor setInversedIndices(SimpleTensor eps) {
        return Tensors.simpleTensor(eps.getName(), eps.getIndices().getInverse());
    }

    //todo static can cause to an indeterminate behavior if CC.resetTensorNames() invoked
    private static final TIntObjectHashMap<Expression[]> cachedSubstitutions = new TIntObjectHashMap<>();

    private static Map<IntArray, Boolean> getEpsilonSymmetries(int indicesSize) {
        Map<IntArray, Boolean> symmetries = cachedLeviCivitaSymmetries.get(indicesSize);
        if (symmetries != null)
            return symmetries;
        symmetries = new HashMap<>();
        Symmetries ss = SymmetriesFactory.createSymmetries(indicesSize);
        ss.addUnsafe(new Symmetry(Combinatorics.createTransposition(indicesSize, 0, 1), true));
        if (indicesSize % 2 == 0)
            ss.addUnsafe(new Symmetry(Combinatorics.createCycle(indicesSize), true));
        else
            ss.addUnsafe(new Symmetry(Combinatorics.createCycle(indicesSize), false));
        for (Symmetry symmetry : ss)
            symmetries.put(symmetry.getPermutation(), symmetry.isAntiSymmetry());
        cachedLeviCivitaSymmetries.put(indicesSize, symmetries);
        return symmetries;
    }

    private static TIntObjectHashMap<Map<IntArray, Boolean>> cachedLeviCivitaSymmetries = new TIntObjectHashMap<>();


}
