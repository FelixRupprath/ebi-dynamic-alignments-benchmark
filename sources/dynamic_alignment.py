from itertools import combinations_with_replacement

from pm4py import ProcessTree
from pm4py.objects.process_tree.obj import Operator


ProcessTree.labels = None


def add_labels_to_process_tree(process_tree: ProcessTree) -> set[str]:
    if process_tree.labels is None:
        if process_tree.operator is None:  # Leaf node
            process_tree.labels = {process_tree.label} if process_tree.label is not None else set()
        else:
            labels = set()
            for child in process_tree.children:
                labels.update(add_labels_to_process_tree(child))
            process_tree.labels = labels
    return process_tree.labels


def _check_unique_labels(process_tree: ProcessTree) -> None:
    """Raise ValueError if any activity label appears in more than one subtree of an operator node.
    τ (silent activity, represented as None) is exempt and may appear multiple times.
    Must be called after add_labels_to_process_tree has populated the .labels cache."""
    if process_tree.operator is None:
        return
    seen: set[str] = set()
    for child in process_tree.children:
        overlap = seen & child.labels
        if overlap:
            raise ValueError(
                f"Process tree does not have unique labels: {overlap!r} occur in multiple subtrees. "
                "dyn_align requires the unique-label property."
            )
        seen |= child.labels
        _check_unique_labels(child)


def check_unique_labels(process_tree: ProcessTree) -> bool:
    if process_tree.operator is None:
        return True
    seen: set[str] = set()
    for child in process_tree.children:
        overlap = seen & child.labels
        if overlap:
            return False
        seen |= child.labels
        if not check_unique_labels(child):
            return False
    return True


def dyn_align(trace: tuple[str, ...], process_tree: ProcessTree) -> int:
    add_labels_to_process_tree(process_tree)  # Initialize label cache
    _check_unique_labels(process_tree)

    cost_table: dict[tuple[tuple[str, ...], int], int] = {}
    # alignment_table = {}

    def _dyn_align(sub_trace: tuple[str, ...], sub_tree: ProcessTree, brute_force: bool = False) -> int:
        if (key := (sub_trace, id(sub_tree))) in cost_table:
            return cost_table[key]

        class Subsequence:
            def __init__(self, subsequence: list[list[str]] | None = None, log_skips: int = 0):
                self.subsequence: list[list[str]] = subsequence or [[]]
                self.log_skips: int = log_skips

            def get_subsequence(self) -> list[list[str]]:
                return [[event for event in sub] for sub in self.subsequence]

            def get_log_skips(self) -> int:
                return self.log_skips

        match sub_tree.operator:
            case None:
                # Base case: leaf node
                cost: int = len(sub_trace) - (1 if sub_tree.label in sub_trace else -1 if sub_tree.label else 0)
                cost_table[key] = cost
                return cost

            case Operator.SEQUENCE:
                if brute_force:
                    # Brute force approach:
                    # Try all possible partitions of the trace and compute the cost of aligning each
                    # The number of partitions is (len(sub_trace) + len(sub_tree.children) - 1) choose len(sub_trace)
                    cost: int = min(sum(_dyn_align(sub_trace[partition[i]:
                                                             partition[i + 1] if i < len(sub_tree.children) - 1
                                                             else len(sub_trace)], child)
                                        for i, child in enumerate(sub_tree.children))
                                    for partition in combinations_with_replacement(range(len(sub_trace) + 1),
                                                                                   len(sub_tree.children) - 1))
                    cost_table[key] = cost
                    return cost

                if len(sub_trace) > 0:
                    # For unique labels in the process tree, we can determine the corresponding tree child for each event
                    # beforehand, and, therefore, len(sub_trace) - floor((len(sub_trace) - 1) / len(sub_tree.children)) is
                    # an upper bound for the number of relevant partitions

                    # Create a dictionary to map events to their corresponding child indices
                    event_to_child_index: dict[str, int] = {event: i
                                                            for i, child in enumerate(sub_tree.children)
                                                            for event in child.labels}

                    sub_traces: dict[int, list[Subsequence]] = {}
                    last_relevant_index: dict[int, int] = {}

                    for i, event in enumerate(sub_trace):
                        if (phi := event_to_child_index.get(event, -1)) >= 0:
                            added = False
                            ref_idx = 0
                            if phi in sub_traces:
                                ref_idx = last_relevant_index[phi]
                                last_relevant_index[phi] = i
                                for subsequence in sub_traces[phi]:
                                    subsequence.subsequence[-1].append(event)
                                    subsequence.log_skips += i - ref_idx - 1
                                added = True
                            for j in reversed(range(phi)):
                                if j in sub_traces and last_relevant_index[j] >= ref_idx:
                                    ref_idx = last_relevant_index[j]
                                    last_relevant_index[phi] = i
                                    sub_traces[phi] = sub_traces.get(phi, []) + [
                                        Subsequence(subsequence.get_subsequence() + [[]] * (phi-j-1) + [[event]],
                                                    subsequence.get_log_skips() + i - ref_idx - 1)
                                        for subsequence in sub_traces[j]]
                                    added = True
                            if not added:
                                sub_traces[phi] = [Subsequence([[]] * phi + [[event]], i)]
                                last_relevant_index[phi] = i

                    ref_idx = 0
                    relevant_subsequences: list[Subsequence] = []
                    for i in reversed(range(len(sub_tree.children))):
                        if i in sub_traces and last_relevant_index[i] >= ref_idx:
                            relevant_subsequences += [
                                Subsequence(subsequence.get_subsequence() + [[]] * (len(sub_tree.children) - i - 1),
                                            subsequence.get_log_skips() + len(sub_trace) - last_relevant_index[i] - 1)
                                for subsequence in sub_traces[i]]
                            ref_idx = last_relevant_index[i]

                    if len(relevant_subsequences) == 0:
                        relevant_subsequences = [Subsequence([[]] * (len(sub_tree.children)), len(sub_trace))]

                    cost: int = min(subsequence.log_skips + sum(_dyn_align(tuple(subsequence.subsequence[i]), child)
                                                                for i, child in enumerate(sub_tree.children))
                                    for subsequence in relevant_subsequences)
                    cost_table[key] = cost
                    return cost
                else:
                    cost: int = sum(_dyn_align(tuple(), child) for child in sub_tree.children)
                    cost_table[key] = cost
                    return cost

            case Operator.XOR:
                # Compute minimum cost of aligning trace with each child
                cost: int = min(_dyn_align(sub_trace, child) for child in sub_tree.children)
                cost_table[key] = cost
                return cost

            case Operator.PARALLEL:
                # Split trace into subtraces for each child and sum up costs of aligning each subtrace with the
                # corresponding child and add cost for unmatched events
                cost: int = len(tuple(event for event in sub_trace if event not in sub_tree.labels))
                for child in sub_tree.children:
                    cost += _dyn_align(tuple(event for event in sub_trace if event in child.labels), child)
                cost_table[key] = cost
                return cost

            case Operator.LOOP:
                do_child = sub_tree.children[0]
                redo_child = sub_tree.children[1]
                n = len(sub_trace)
                loop_cost_table: dict[tuple[int, int, int], int] = {}

                # Precompute the fixed cost of executing the redo-part on an empty sub-trace.
                cost_redo_tau: int = _dyn_align((), redo_child)

                def _loop_cost(start: int, end: int, child: ProcessTree) -> int:
                    loop_key = (start, end, id(child))
                    if loop_key in loop_cost_table:
                        return loop_cost_table[loop_key]
                    result = _dyn_align(sub_trace[start:end], child)
                    loop_cost_table[loop_key] = result
                    return result

                # Compute trivial upper bound by aligning the whole sub_trace to the do_child.
                # Leave the computation early if a perfect alignment was achieved.
                if (upper_bound := _loop_cost(0, n, do_child)) == 0:
                    cost_table[key] = 0
                    return 0

                # Iteratively compute the minimal cost at each position in the sub_trace, but only consider
                # intermediate results that can actually improve against the upper bound.
                # Finally, return the minimal cost at the last position.
                intermediate_cost_after_do: list[int] = []
                intermediate_cost_after_redo: list[int] = [0]

                for i in range(n + 1):
                    best_do = upper_bound + 1
                    for j in range(i + 1):
                        if intermediate_cost_after_redo[j] < best_do:
                            if (cur_do := intermediate_cost_after_redo[j] + _loop_cost(j, i, do_child)) < best_do:
                                if (best_do := cur_do) == 0:
                                    break
                    intermediate_cost_after_do.append(best_do)

                    if i == n:
                        break

                    intermediate_cost_after_redo[-1] = min(intermediate_cost_after_redo[-1],
                                                           best_do + cost_redo_tau,
                                                           )

                    best_redo = upper_bound + 1
                    for j in range(i + 1):
                        if intermediate_cost_after_do[j] < best_redo:
                            if (cur_redo := intermediate_cost_after_do[j] + _loop_cost(j, i + 1, redo_child)) < best_redo:
                                if (best_redo := cur_redo) == 0:
                                    break
                    intermediate_cost_after_redo.append(best_redo)

                cost: int = intermediate_cost_after_do[n]
                cost_table[key] = cost
                return cost

            case Operator.OR:
                cost: int = len(sub_trace)
                gain: list[int] = []
                for child in sub_tree.children:
                    split_trace: tuple[str, ...] = tuple(event for event in sub_trace if event in child.labels)
                    gain.append(len(split_trace) - _dyn_align(split_trace, child))
                if max(gain) > 0:
                    for g in gain:
                        if g > 0:
                            cost -= g
                else:
                    cost -= max(gain)
                cost_table[key] = cost
                return cost

            case _:
                raise NotImplementedError(f"Operator {sub_tree.operator} not supported!")

    return _dyn_align(trace, process_tree)
