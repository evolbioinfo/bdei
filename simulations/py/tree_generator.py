import logging
import random

import numpy as np
from ete3 import TreeNode

TIME_TILL_NOW = 'T'
DIST_TO_START = 'D'
STATE = 'state'


def simulate_tree_gillespie(model, max_time=np.inf, max_sampled=np.inf,
                            state_feature=STATE, root_state=None, state_frequencies=None):
    """
    Simulates the tree evolution from a root over the given time based on the given model.

    :param root_state: root state, if not specified,
        will be drawn randomly from model states according to their equilibrium frequencies.
    :param state_frequencies: array of equilibrium frequencies of model states
        (to be used to draw the root state if it is not given). By default all equal.
    :param max_sampled: maximal number of sampling node (when reached, the simulation stops)
    :param max_time: float, time over which we generate a tree.
    :param model: bd_models.Model
    :return: the simulated tree (ete3.Tree).
    """
    num_states = len(model.states)
    if root_state is None:
        if state_frequencies is None:
            state_frequencies = [1 / num_states] * num_states
        root_state = np.random.choice(model.states, size=1, p=state_frequencies)[0]
    # evolve till the time is up, following Gillespie
    time = 0
    infectious_nums = np.zeros(num_states, dtype=np.int)
    infectious_nums[root_state.index] = 1
    sampled_nums = np.zeros(num_states, dtype=np.int)

    infectious_state2id = [set() for _ in model.states]
    cur_id = 0, 0
    infectious_state2id[root_state.index].add(cur_id)
    id2time = {}
    id2parent = {}
    sampled_id2state = {}

    while infectious_nums.sum() and sampled_nums.sum() < max_sampled and time < max_time:
        # first we need to calculate rate sum
        rate_sums = model.rates.dot(infectious_nums)
        total_rate = rate_sums.sum()

        # now let us see when next event takes place
        time += np.random.exponential(1 / total_rate, 1)[0]

        # Check if the time is up
        if time > max_time:
            break

        # now let us see which event will happen
        random_event = np.random.uniform(0, 1, 1)[0] * total_rate

        # case 1: state transition
        if random_event < rate_sums[0]:
            transition_rates = model.rates[0, :] * infectious_nums
            for i in range(num_states):
                if random_event < transition_rates[i]:
                    state = model.states[i]
                    infectious_nums[state.index] -= 1
                    infectious_nums[state.next_state.index] += 1

                    infectious_state2id[state.next_state.index].add(random_pop(infectious_state2id[state.index]))
                    break
                random_event -= transition_rates[i]
            continue
        random_event -= rate_sums[0]

        # case 2: transmission
        if random_event < rate_sums[1]:
            transmission_rates = model.rates[1, :] * infectious_nums
            for i in range(num_states):
                if random_event < transmission_rates[i]:
                    state = model.states[i]
                    infectious_nums[state.recipient.index] += 1

                    cur_id = cur_id[0] + 1, 0
                    parent_id = random_pop(infectious_state2id[state.index])
                    donor_id = parent_id[0], parent_id[1] + 1
                    infectious_state2id[state.index].add(donor_id)
                    infectious_state2id[state.recipient.index].add(cur_id)
                    id2parent[cur_id] = parent_id
                    id2parent[donor_id] = parent_id
                    id2time[parent_id] = time
                    break
                random_event -= transmission_rates[i]
            continue
        random_event -= rate_sums[1]

        # case 3: sampling
        sampling_rates = model.rates[2, :] * infectious_nums
        for i in range(num_states):
            if random_event < sampling_rates[i]:
                state = model.states[i]
                infectious_nums[state.index] -= 1

                sampled_id = random_pop(infectious_state2id[state.index])
                id2time[sampled_id] = time

                if np.random.uniform(0, 1, 1)[0] < model.ps[state.index]:
                    sampled_id2state[sampled_id] = state
                    sampled_nums[state.index] += 1

                break
            random_event -= sampling_rates[i]
    if max_time == np.inf:
        max_time = time

    return reconstruct_tree(id2parent, id2time, sampled_id2state, max_time, state_feature=state_feature)


def reconstruct_tree(id2parent, id2time, sampled_id2state, max_time, state_feature=STATE):
    if not sampled_id2state:
        return None

    root = None
    id2node = {}
    for id, state in sampled_id2state.items():
        time = id2time[id]
        node = TreeNode(dist=time - (0 if not id in id2parent else id2time[id2parent[id]]), name=id[0])
        node.add_feature(DIST_TO_START, time)
        node.add_feature(state_feature, state)
        id2node[id] = node
        while node is not None and id in id2parent:
            parent_id = id2parent[id]
            if parent_id in id2node:
                id2node[parent_id].add_child(node)
                break
            parent_time = id2time[parent_id]
            parent = TreeNode(dist=parent_time - (0 if not parent_id in id2parent else id2time[id2parent[parent_id]]))
            parent.add_feature(DIST_TO_START, parent_time)
            id2node[parent_id] = parent
            parent.add_child(node)
            node, id = parent, parent_id
        if id not in id2parent:
            root = node
    # remove internal nodes with just one child
    for node in root.traverse('postorder'):
        if len(node.children) == 1:
            child = node.children[0]
            child.dist += node.dist
            if not node.is_root():
                parent = node.up
                parent.remove_child(node)
                parent.add_child(child)
            else:
                root = child
                child.up = None
    # annotate time till now
    for node in root.traverse('postorder'):
        node.add_feature(TIME_TILL_NOW, max_time - getattr(node, DIST_TO_START))
    return root


def random_pop(elements):
    """
    Removes a random element from a list and returns it.
    :param elements: list of elemetns
    :return: the selected element
    """
    element = random.sample(elements, 1)[0]
    elements.remove(element)
    return element


def generate_forest(model, max_time=np.inf, min_tips=1_000, max_sampled=np.inf, keep_nones=False, state_feature=STATE,
                    root_state=None, state_frequencies=None):
    total_n_tips = 0
    forest = []
    total_trees = 0
    sampled_trees = 0
    while total_n_tips < min_tips:
        tree = simulate_tree_gillespie(model, max_time=max_time, max_sampled=max_sampled,
                                       state_feature=state_feature,
                                       root_state=root_state, state_frequencies=state_frequencies)
        total_trees += 1
        if tree:
            total_n_tips += len(tree)
            sampled_trees += 1
        if tree or keep_nones:
            forest.append(tree)

    logging.info('Total number of tips n={}.'.format(total_n_tips))
    return forest
