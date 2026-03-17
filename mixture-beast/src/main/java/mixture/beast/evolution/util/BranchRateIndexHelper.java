package mixture.beast.evolution.util;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.parameter.RealParameter;

import java.util.Arrays;

/**
 * Shared helper that builds the deterministic mapping:
 *
 *   branch (identified by child node nr) -> index in the shared rawRates vector
 *
 * using the exact same convention everywhere:
 *   - root child index = -1
 *   - all non-root branches are assigned indices 0..(nodeCount-2)
 *     in ascending child-node nr order
 *
 * This matches the deterministic mapping used in RelaxedRatesPriorSVS.
 */
public final class BranchRateIndexHelper {

    private BranchRateIndexHelper() {
        // utility class
    }

    public static final class Mapping {
        private final int nodeCount;
        private final int rootNr;
        private final int[] idxMap;

        private Mapping(final int nodeCount, final int rootNr, final int[] idxMap) {
            this.nodeCount = nodeCount;
            this.rootNr = rootNr;
            this.idxMap = idxMap;
        }

        public int getNodeCount() {
            return nodeCount;
        }

        public int getRootNr() {
            return rootNr;
        }

        public boolean matches(final Tree tree) {
            return tree != null
                    && tree.getRoot() != null
                    && tree.getNodeCount() == nodeCount
                    && tree.getRoot().getNr() == rootNr;
        }

        public int idxForNodeNr(final int nr) {
            if (nr < 0 || nr >= nodeCount) {
                throw new IllegalArgumentException("Node nr out of range: " + nr + " (nodeCount=" + nodeCount + ")");
            }
            return idxMap[nr];
        }

        public int idxForNode(final Node node) {
            if (node == null) {
                throw new IllegalArgumentException("Node is null");
            }
            return idxForNodeNr(node.getNr());
        }

        public int[] copyIndexMap() {
            return Arrays.copyOf(idxMap, idxMap.length);
        }
    }

    public static Mapping buildDeterministic(final Tree tree) {
        if (tree == null) {
            throw new IllegalArgumentException("Tree is null");
        }
        if (tree.getRoot() == null) {
            throw new IllegalArgumentException("Tree root is null");
        }

        final int nNodes = tree.getNodeCount();
        final int rootNr = tree.getRoot().getNr();

        if (nNodes < 1) {
            throw new IllegalArgumentException("Tree must contain at least one node");
        }
        if (rootNr < 0 || rootNr >= nNodes) {
            throw new IllegalArgumentException("Root nr out of range: " + rootNr + " (nodeCount=" + nNodes + ")");
        }

        final Node[] byNr = new Node[nNodes];
        for (int i = 0; i < nNodes; i++) {
            final Node node = tree.getNode(i);
            if (node == null) {
                throw new IllegalStateException("tree.getNode(" + i + ") returned null");
            }

            final int nr = node.getNr();
            if (nr < 0 || nr >= nNodes) {
                throw new IllegalArgumentException("Node nr out of range: " + nr + " (nodeCount=" + nNodes + ")");
            }
            if (byNr[nr] != null) {
                throw new IllegalStateException("Duplicate node nr encountered: " + nr);
            }
            byNr[nr] = node;
        }

        for (int nr = 0; nr < nNodes; nr++) {
            if (byNr[nr] == null) {
                throw new IllegalStateException("Missing node for nr=" + nr);
            }
        }

        final int[] idxMap = new int[nNodes];
        Arrays.fill(idxMap, -2);

        int k = 0;
        for (int nr = 0; nr < nNodes; nr++) {
            if (nr == rootNr) {
                idxMap[nr] = -1;
            } else {
                idxMap[nr] = k++;
            }
        }

        if (k != nNodes - 1) {
            throw new IllegalStateException("Internal error: expected to map " + (nNodes - 1)
                    + " non-root nodes but mapped " + k);
        }

        return new Mapping(nNodes, rootNr, idxMap);
    }

    public static void validateRatesDimension(final Tree tree,
                                              final RealParameter rates,
                                              final String ownerName) {
        if (tree == null) {
            throw new IllegalArgumentException(ownerName + ": tree is null");
        }
        if (rates == null) {
            throw new IllegalArgumentException(ownerName + ": rates is null");
        }

        final int expected = tree.getNodeCount() - 1;
        if (rates.getDimension() != expected) {
            throw new IllegalArgumentException(ownerName + ": rates must have dimension (nodeCount - 1). Found "
                    + rates.getDimension() + " vs " + expected);
        }
    }

    public static Mapping ensureUpToDate(final Tree tree,
                                         final Mapping mapping,
                                         final RealParameter rates,
                                         final String ownerName) {
        validateRatesDimension(tree, rates, ownerName);
        if (mapping == null || !mapping.matches(tree)) {
            return buildDeterministic(tree);
        }
        return mapping;
    }
}
