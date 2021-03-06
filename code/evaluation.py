#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by sunder on 2017/12/12

import json
import logging

import numpy as np

import data_prepare

logger = logging.getLogger('mylogger')


def calc_advantage2(actions, standard_ys, config):
    rewards = []
    numbers = []
    gold_numbers = []
    for action, y in zip(actions, standard_ys):
        r = []
        n = 0
        visited = set()
        gold_triples = _triplelist2triples_(y, config, strict=True)
        gold_numbers.append(len(gold_triples))
        triple_list = list(action)
        predict_triples = [tuple(triple_list[i:i + 3]) for i in range(0, len(triple_list), 3)]
        for i, triple in enumerate(predict_triples):
            if i < len(gold_triples):
                if triple in gold_triples:
                    if triple in visited:
                        r.append(0.0)
                        r.append(0.0)
                        r.append(0.0)
                    else:
                        visited.add(triple)
                        r.append(1.0)
                        r.append(1.0)
                        r.append(1.0)
                        n += 1
                else:
                    r.append(0.0)
                    r.append(0.0)
                    r.append(0.0)
            elif triple == config.NA_TRIPLE:
                r.append(0.5)
                r.append(0.5)
                r.append(0.5)
            else:
                r.append(0.0)
                r.append(0.0)
                r.append(0.0)
        rewards.append(r)
        numbers.append(n)
    logger.debug(map(str, actions[0]))
    logger.debug(map(str, standard_ys[0]))
    logger.debug(map(str, rewards[0]))
    logger.debug('%d / %d' % (numbers[0], gold_numbers[0]))
    return rewards, numbers



def calc_advantage(actions, standard_ys, config):
    rewards = []
    numbers = []
    for action, y in zip(actions, standard_ys):
        # use the f1 as reward
        f1, precision, recall, correct_num, predict_number, gold_number = compare_([action], [y], name='rewards',
                                                                                   config=config, is_show=False)
        reward = f1
        rewards.append(reward)
        numbers.append(predict_number)
    return rewards, numbers


def compare(predict, gold, config,  result_file,show_rate=None, simple=True):
    normal_triples_gold = []  # normal triples
    normal_triples_predict = []  # normal triples
    multi_label_gold = []  # multi label triples
    multi_label_predict = []  # multi label triples
    over_lapping_gold = []  # overlapping triples
    over_lapping_predict = []  # overlapping triples
    is_relation_first = True
    for p, g in zip(predict, gold):
        if data_prepare.is_normal_triple(g, is_relation_first):
            normal_triples_gold.append(g)
            normal_triples_predict.append(p)
        if data_prepare.is_multi_label(g, is_relation_first):
            multi_label_gold.append(g)
            multi_label_predict.append(p)
        if data_prepare.is_over_lapping(g, is_relation_first):
            over_lapping_gold.append(g)
            over_lapping_predict.append(p)
    f1, precision, recall, _, _, _ = compare_( predict, gold, 'ALL', config,  result_file, show_rate, strict=True)
    #if simple:
        #return f1, precision, recall
    compare_(normal_triples_predict, normal_triples_gold, 'Normal-Triples', config, show_rate)
    compare_(multi_label_predict, multi_label_gold, 'Multi-Label', config, show_rate)
    compare_(over_lapping_predict, over_lapping_gold, 'Over-Lapping', config, show_rate)

    # sentences contains 1, 2, 3, 4, and >5 triples
    triples_size_1_gold, triples_size_2_gold, triples_size_3_gold, triples_size_4_gold, triples_size_5_gold = [], [], [], [], []
    triples_size_1_predict, triples_size_2_predict, triples_size_3_predict, triples_size_4_predict, triples_size_5_predict = [], [], [], [], []
    for p, g in zip(predict, gold):
        g_triples = set([tuple(g[i:i + 3]) for i in range(0, len(g), 3)])
        if len(g_triples) == 1:
            triples_size_1_predict.append(p)
            triples_size_1_gold.append(g)
        elif len(g_triples) == 2:
            triples_size_2_predict.append(p)
            triples_size_2_gold.append(g)
        elif len(g_triples) == 3:
            triples_size_3_predict.append(p)
            triples_size_3_gold.append(g)
        elif len(g_triples) == 4:
            triples_size_4_predict.append(p)
            triples_size_4_gold.append(g)
        else:
            triples_size_5_predict.append(p)
            triples_size_5_gold.append(g)
    compare_(triples_size_1_predict, triples_size_1_gold, 'Sentence-1-Triple', config, show_rate)
    compare_(triples_size_2_predict, triples_size_2_gold, 'Sentence-2-Triple', config, show_rate)
    compare_(triples_size_3_predict, triples_size_3_gold, 'Sentence-3-Triple', config, show_rate)
    compare_(triples_size_4_predict, triples_size_4_gold, 'Sentence-4-Triple', config, show_rate)
    compare_(triples_size_5_predict, triples_size_5_gold, 'Sentence-5-Triple', config, show_rate)
    return f1, precision, recall


def _triplelist2triples_(triple_list, config, strict=True, keep_order=False):
    """
     >>> _triplelist2triples_([1,2,3, 2,5,0])
     {(1,2,3),(2,5,0)}
     >>> _triplelist2triples_([1,2,3, 1,2,3, 2,5,0])
     {(1,2,3),(2,5,0)}
     >>> _triplelist2triples_([1,2,3, 2,5,0].extend(config.NA_TRIPLE))
     {(1,2,3),(2,5,0)}
    """
    triple_list = list(triple_list)
    if keep_order:
        triples = [tuple(triple_list[i:i + 3]) for i in range(0, len(triple_list), 3)]
    else:
        triples = set([tuple(triple_list[i:i + 3]) for i in range(0, len(triple_list), 3)])
    if strict:
        # Remove the triple exactly the same with NA triple
        while config.NA_TRIPLE in triples:
            triples.remove(config.NA_TRIPLE)
    else:
        # The triple with NA relation is regraded as NA triple
        valid_triples = set()
        for triple in triples:
            if triple[0] != config.NA_TRIPLE[0]:
                valid_triples.add(triple)
        triples = valid_triples

    return triples


def triples2entities(triples):
    """
    :param triples:
    :return:
    >>> triples2entities([[1,2,3], [0, 3,4]])
    [2,3,4]
    >>> triples2entities([[1,2,3], [1,2,3]])
    [2,3]
    """
    entities = []
    for triple in triples:
        entities.extend(triple[1:])
    return list(set(entities))


def triples2relations(triples):
    """
    :param triples:
    :return:
    >>> triples2relations([[1,2,3], [0, 3,4]])
    [1,0]
    >>> triples2relations([[1,2,3], [1,2,3]])
    [1]
    """
    relations = []
    for triple in triples:
        relations.append(triple[0])
    return list(set(relations))


def error_analyse(predicts, gold, config, entity_or_relation='entity'):
    predict_number = 0
    gold_number = 0
    correct_num = 0
    func = triples2entities if entity_or_relation == 'entity' else triples2relations
    for p, g in zip(predicts, gold):
        p_triples = _triplelist2triples_(p, config)
        g_triples = _triplelist2triples_(g, config)
        p_elements = func(p_triples)
        g_elements = func(g_triples)
        predict_number += len(p_elements)
        gold_number += len(g_elements)
        result = [1 if e in g_elements else 0 for e in p_elements]
        correct_num += sum(result)

    logger.debug('Error Analyse: %s: Predict number %d, Gold number %d, Correct number %d' % (
    entity_or_relation, predict_number, gold_number, correct_num))
    precision = correct_num * 1.0 / predict_number if predict_number > 0 else 0.
    recall = correct_num * 1.0 / gold_number if gold_number > 0 else 0.0
    f1 = 2 * precision * recall / (precision + recall) if precision * recall > 0 else 0.
    logger.info('Error Analyse: %s: Precision %s, Recall %s, F1 %s' % (entity_or_relation, precision, recall, f1))


def compare_(predict, gold, name, config, result_file, show_rate=None, is_show=True, strict=True):
    predict_number = 0
    gold_number = 0
    correct_num = 0
    for p, g in zip(predict, gold):
        p_triples = _triplelist2triples_(p, config, strict=strict)
        g_triples = _triplelist2triples_(g, config, strict=strict)
        predict_number += len(p_triples)
        gold_number += len(g_triples)

        result = [1 if p_t in g_triples else 0 for p_t in p_triples]
        correct_num += sum(result)

        logging.basicConfig(level=logging.DEBUG,  # 控制台打印的日志级别
                            filename='new.log',
                            filemode='a',  ##模式，有w和a，w就是写模式，每次都会重新写日志，覆盖之前的日志
                            # a是追加模式，默认如果不写的话，就是追加模式
                            format=
                            '%(asctime)s - %(pathname)s[line:%(lineno)d] - %(levelname)s: %(message)s'
                            # 日志格式
                            )



        if np.random.uniform() < show_rate:
            logger.debug('%s: predict %s' % (name, p))
            logger.debug('%s: gold    %s' % (name, g))
            logger.debug(
                '%s: ----------------------------------------------------------------------------- result %s/%s' % (
                name, sum(result), len(g_triples)))
    # write deocde result to files !
    np.savez(predict, result_file)
    precision = correct_num * 1.0 / predict_number if predict_number > 0 else 0.
    recall = correct_num * 1.0 / gold_number if gold_number > 0 else 0.0
    f1 = 2 * precision * recall / (precision + recall) if precision * recall > 0 else 0.
    if is_show:
        logger.info('%s: Instance number %d' % (name, len(gold)))
        logger.info(
            '%s: Predict number %d, Gold number %d, Correct number %d' % (
            name, predict_number, gold_number, correct_num))
        logger.info('%s: Precision %.3f, Recall %.3f, F1 %.3f' % (name, precision, recall, f1))
    return f1, precision, recall, correct_num, predict_number, gold_number


def sent_id2sent_str(sent_id, id2words):
    words = []
    for idx in sent_id:
        if idx > 0:
            try:
                word = id2words[idx]
            except:
                word = 'None-%d' % idx
            words.append(word)
    sent_str = ' '.join(words).encode('utf-8').strip()
    return sent_str


def triple_id2triple_str(triple_id, sent_id, id2words, id2relations, is_relation_first, config):
    assert len(triple_id) == 3
    entity_1_str, entity_2_str, relation_str = 'None', 'None', 'None'
    if is_relation_first:
        r_id, e_1_position_id, e_2_position_id = triple_id[0], triple_id[1], triple_id[2]
    else:
        r_id, e_1_position_id, e_2_position_id = triple_id[2], triple_id[0], triple_id[1]
    if e_1_position_id < config.max_sentence_length:
        try:
            entity_1_str = id2words[sent_id[e_1_position_id]]
        except:
            entity_1_str = 'None-%s-%s' % (e_1_position_id, sent_id[e_1_position_id])
    if e_2_position_id < config.max_sentence_length:
        try:
            entity_2_str = id2words[sent_id[e_2_position_id]]
        except:
            entity_2_str = 'None-%s-%s' % (e_1_position_id, sent_id[e_1_position_id])
    if r_id < config.relation_number:
        try:
            relation_str = id2relations[r_id]
        except:
            relation_str = 'None-%s' % r_id
    return '[%s, %s, %s]' % (
    entity_1_str.encode('utf-8').strip(), entity_2_str.encode('utf-8').strip(), relation_str.encode('utf-8').strip())


def triples2triples_str(triples, sent_id, id2words, id2relations, is_relation_first, config):
    triples_str = []
    for triple in triples:
        triple_string = triple_id2triple_str(triple, sent_id, id2words, id2relations, is_relation_first, config)
        triples_str.append(triple_string)
    return '\t'.join(triples_str)


def _reverse_dict_(a_dict):
    new_dict = {v: k for k, v in a_dict.items()}
    return new_dict


def visualize(sents_id, gold, predict, files_name, config, is_relation_first=True):
    print ('Visualizing ...')
    print (config.words2id_filename)
    print (config.relations2id_filename)
    words2id = json.load(open(config.words2id_filename, 'r'))
    relations2id = json.load(open(config.relations2id_filename, 'r'))
    id2words = _reverse_dict_(words2id)
    id2relations = _reverse_dict_(relations2id)

    f1 = open(files_name[0], 'w')
    f2 = open(files_name[1], 'w')
    f3 = open(files_name[2], 'w')
    for d, g, p in zip(sents_id, gold, predict):
        if data_prepare.is_normal_triple(g, is_relation_first):
            f = f1
        elif data_prepare.is_multi_label(g, is_relation_first):
            f = f2
        else:
            f = f3

        f.write(sent_id2sent_str(d, id2words))
        f.write('\n')
        g_triples = _triplelist2triples_(g, config, keep_order=True)
        p_triples = _triplelist2triples_(p, config, keep_order=True)
        g_triples_string = triples2triples_str(g_triples, d, id2words, id2relations, is_relation_first, config)
        p_triples_string = triples2triples_str(p_triples, d, id2words, id2relations, is_relation_first, config)
        f.write('Gold:   \t' + g_triples_string)
        f.write('\n')
        f.write('Predict:\t' + p_triples_string)
        f.write('\n\n')
    f1.close()
    f2.close()
    f3.close()


def static_every_time_step(batch_actions, predict_space_size, copy_space_size, max_step):
    records = []
    for step in range(max_step):
        if step % 3 == 0:
            record = dict([(i, 0) for i in range(predict_space_size)])
        else:
            record = dict([(i, 0) for i in range(copy_space_size)])
        records.append(record)
    for actions in batch_actions:
        for step, act in enumerate(actions):
            if step < max_step:
                count = records[step].get(act, 0)
                records[step][act] = count + 1
    records = [sorted(record.items(), key=lambda x: x[0]) for record in records]
    records = [[x[1] for x in record] for record in records]
    # records = [[x if x < 50 else 50 for x in record] for record in records]
    records = [[x * 1.0 / sum(record) for x in record] for record in records]
    return records





if __name__ == '__main__':
    import doctest

    doctest.testmod()
