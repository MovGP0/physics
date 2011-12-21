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
package redberryphysics.core.oneloop;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import redberry.core.context.CC;
import redberry.core.tensor.Expression;
import redberry.core.tensor.SimpleTensor;
import redberry.core.tensor.Tensor;
import redberry.core.tensor.TensorIterator;
import redberry.core.transformation.ExpandBrackets;
import redberry.core.transformation.RenameConflictingIndexes;
import redberry.core.transformation.Transformer;
import redberry.core.transformation.collect.CollectTerms;
import redberry.core.transformation.collect.PatternSplitCriteria;
import redberry.core.transformation.collect.SplitPattern;
import redberry.core.transformation.contractions.IndexesContractionsTransformation;

/**
 *
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
public class R2Gravity {
    public static final Expression Ric =
            new Expression("Ric_ab=E^r_a*E^d_c*R^c_bdr");
    public static final Expression Riman =
            new Expression("R^a_bcd=D[w^a_db,x^c]-D[w^a_cb,x^d]+w^a_cr*w^r_db-w^a_dr*w^r_cb");
    public static final Expression Torsion =
            new Expression("T^a_bc=D[h^a_c,x^b]-D[h^a_b,x^c]+w^a_bd*e^d_c-w^a_cd*e^d_b");
    public static final Expression eTetrad =
            new Expression("e^a_b=d^a_b+h^a_b");
    public static final Expression ETetrad =
            new Expression("E^a_b=d^a_b-h^a_b+h^a_c*h^c_b");
    public static final Expression metricUP =
            new Expression("G^ab=g^ab-g^ca*h^b_c-g^cb*h^a_c+g^cb*h^a_d*h^d_c+g^ca*h^b_d*h^d_c+g^cd*h^a_c*h^b_d");
    public static final Expression sqrt =
            new Expression("sqrt=1+h^a_a+(1/2)*(h^s_s*h^l_l-h^s_l*h^l_s)");
    public static final Expression tetradGaugeFix =
            new Expression("Gf_a=f1*h^b_a*p_b+f2*g^pq*g_ab*h^b_q*p_p+f3*h^q_q*p_a");
    public static final Expression Lagrangian =
            new Expression("Lagrange = sqrt*(g^ab*Ric_ab+(e1*g_ab*G^xp*G^yq+e2*E^x_a*E^p_b*G^yq+e3*E^x_b*E^p_a*G^yq)*T^a_xy*T^b_pq+e6*Ric_ab*Ric_cd*g^ab*g^cd+e5*Ric_ab*Ric_cd*g^ac*g^bd+Gf_a*Gf_b*g^ab)+f*g^pq*g_ab*i*h^b_q*p_p*g^cd*I*h^a_d*p_c");
    public final Expression Lagrange;
    public static final Expression[] substitutionsQueue = {
        Ric, Riman, Torsion, eTetrad, ETetrad, metricUP, sqrt, tetradGaugeFix};
    public static final SimpleTensor h = (SimpleTensor) CC.parse("h^a_b");
    public static final SimpleTensor w = (SimpleTensor) CC.parse("w^a_bc");

    public R2Gravity() {
        Lagrange = Lagrangian.clone();
        Lagrange.eval(substitutionsQueue);
        Lagrange.eval(
                new Transformer(RenameConflictingIndexes.INSTANCE),
                IndexesContractionsTransformation.CONTRACTIONS_WITH_KRONECKER,
                new Transformer(ExpandBrackets.EXPAND_EXCEPT_SYMBOLS),
                IndexesContractionsTransformation.CONTRACTIONS_WITH_KRONECKER);
        int matches;
        TensorIterator it = Lagrange.right().iterator();
        Tensor summand;
        while (it.hasNext()) {
            summand = it.next();
            matches = 0;
            for (Tensor t : summand)
                if (t instanceof SimpleTensor) {
                    int name;
                    if ((name = ((SimpleTensor) t).getName()) == h.getName() || name == w.getName())
                        matches++;
                }
            //leaving only no greate than qubic terms
            if (matches > 3)
                it.remove();
        }

        List<SimpleTensor> collectingMatchers = new ArrayList<>();
        collectingMatchers.add(h);
        collectingMatchers.add(w);

        SplitPattern splitPattern = new SplitPattern(collectingMatchers, Collections.EMPTY_LIST);
        //449 summands
        Lagrange.eval(new CollectTerms(new PatternSplitCriteria(splitPattern)));
    }
}
