import re,json,nltk,random,os,csv,requests
from collections import Counter,defaultdict

# home="E:/text_mining/0Data_Collection/corpus/ncRNA-disease/multi_re3/"

home="E:/text_mining/0Data_Collection/211/decode2/seq2seq/gene/"

# f_path="E:/text_mining/0Data_Collection/ncrna_diease_pair/MNDR/ncRNA_disease"
f_path="E:/text_mining/0Data_Collection/ncrna_diease_pair/train.txt"
# train_filtered_filename="E:/text_mining/0Data_Collection/ncrna_diease_pair/train_filtered.json"

gene_disease=home+"train_data01"
# train_filtered_filename=home+"train_filtered.json"
ncrna_disease=home+"train_data_modified"

test_file=home+"decode"

train_filtered_filename= home+"train_filtered.json"
# train_filtered_filename2= home2+"train_filtered.json"


origin_file_path = os.path.join(home, 'origin/')
origin_train_filename = os.path.join(origin_file_path, 'train.json')
origin_test_filename = os.path.join(origin_file_path, 'test.json')
origin_valid_filename = os.path.join(origin_file_path, 'valid.json')
origin_train_filtered_filename = os.path.join(origin_file_path, 'train_filtered.json')
raw_train_filename = os.path.join(origin_file_path, 'raw_train.json')
raw_test_filename = os.path.join(origin_file_path, 'raw_test.json')
raw_valid_filename = os.path.join(origin_file_path, 'raw_valid.json')
words2id_filename = os.path.join(home, 'words2id.json')
relations2id_filename = os.path.join(home, 'relations2id.json')
relation2count_filename = os.path.join(home, 'relation2count.json')
vector_bin_filename = os.path.join(home, 'nyt_vec.bin')
words_id2vector_filename = os.path.join(home, 'words_id2vector.json')
MAX_SENTENCE_LENGTH = 100

#biobert training data
# biobert_home='E:/text_mining/0Data_Collection/corpus/ncRNA-disease/biobert'
biobert_home='E:/text_mining/0Data_Collection/corpus/gene-disease/biobert'
train_f=os.path.join(biobert_home,"train.tsv")
dev_f=os.path.join(biobert_home,"dev.tsv")
test_f=os.path.join(biobert_home,"test.tsv")
test_f_multi=home+"test.tsv"

def prepress_train_data(filename):
    with open(filename, 'r',encoding="utf-8") as fr,open(home+"train_data_modified", 'w',encoding="utf-8") as fw:
        all,error=0,0

        for line in fr.readlines():
            change_e1,change_e2 = False,False
            e1,e2,text,rel=line.strip().split("\t")
            e1 = e1.strip(",")
            e2 = e2.strip(",")
            e1=e1.strip(":")
            e2=e2.strip(":")
            all+=1
            sent_words = nltk.word_tokenize(text)
            towrite = True
            if e1.split()[-1] not in sent_words :
                e1=e1.split()[-1]
                for word in sent_words:
                    if word.startswith(e1) or word.endswith(e1) :
                        e1=word
                        change_e1=True
                        break
                if not change_e1:
                    e1=e1.split(":")[0]
                    if e1 not in sent_words:
                        for i,word in enumerate(sent_words):
                            if e1 in word and "+" in word:
                                word_s=word.split("+")
                                if e1 in word_s:
                                    sent_words=sent_words[:i]+word_s+sent_words[i+1:]
                                break
                            if e1 in word and "/" in word:
                                word_s=word.split("/")
                                if e1 in word_s:
                                    sent_words=sent_words[:i]+word_s+sent_words[i+1:]
                                break
                        if e1 not in sent_words:
                            error+=1
                            towrite = False
            if e2.split()[-1] not in sent_words:
                e2 = e2.split()[-1]
                for word in sent_words:
                    if word.startswith(e2) or  word.endswith(e2):
                        e2 = word
                        change_e2=True
                        break
                if not change_e2:
                    e2 = e2.split(":")[0]
                    if e2 not in sent_words:
                        for i,word in enumerate(sent_words):
                            if e2 in word and "+" in word:
                                word_s=word.split("+")
                                if e2 in word_s:
                                    sent_words=sent_words[:i]+word_s+sent_words[i+1:]
                                break
                            if e2 in word and "/" in word:
                                word_s=word.split("/")
                                if e2 in word_s:
                                    sent_words=sent_words[:i]+word_s+sent_words[i+1:]
                                break
                        if e2 not in sent_words:
                            error += 1
                            towrite = False
            if towrite:
                fw.write(e1 + '\t' + e2 + '\t' + text + '\t' + rel + '\n')

        print(all,error)
# prepress_train_data(ncrna_disease)


def read_json(filename):
    data = []
    with open(filename, 'r',encoding='utf-8') as f:
        for line in f:
            a_data = json.loads(line)
            data.append(a_data)
    return data

def txt2json(filename):
    types = []
    txt2rel=defaultdict(list)
    same_entity=0
    with open(filename, 'r',encoding="utf-8") as f:
        for line in f.readlines():
            e1,e2,text,rel=line.strip().split("\t")
            # if rel=="SusceptibilityMutation" or rel=="ModifyingMutation" or rel=="Therapeutic":
            #     continue
            # rel=rel.split("[")[0].strip()
            # rel=rel.replace("N/A","None")
            # rel=rel.replace("expressed","Expression")
            # rel = rel.replace("mutated", "Mutation")
            # rel = rel.replace("regulated", "Regulation")
            types.append(rel)
            if text.startswith('"'):
                text=text[1:-1]
            txt2rel[text].append((e1,e2,rel))
            # if (e1,e2,rel) not in txt2rel[text.strip()]:
            #     txt2rel[text].append((e1,e2,rel))
        # print(Counter(types)) Counter({'up-Regulation': 4719, 'down-Regulation': 3030, 'Regulation': 1580, 'Expression': 860, 'differentially Expression': 521, 'None': 341, 'Mutation': 92, 'Interaction': 88, 'Locus': 50, 'Epigenetics': 28, 'Association': 17, 'associated': 8, 'detectable': 1})
        print(Counter(types),len(txt2rel))

        data=[]
        for text,triples in txt2rel.items():
             a_data={}
             a_data["sentText"] = text
             a_data["relationMentions"]=[]
             for triple_ in triples:
                triple={}
                triple["em1Text"] = triple_[0].strip(".")
                triple["em2Text"] = triple_[1].strip(".")
                if triple["em1Text"]==triple["em2Text"]:
                    # print(triple)
                    same_entity+=1
                    continue
                triple["label"] = triple_[2]
                if triple["em1Text"] and triple["em2Text"] and triple not in a_data["relationMentions"]:
                    a_data["relationMentions"].append(triple)
             data.append(a_data)
        print("errorous examples %d"%same_entity )
        return data


#   flag is used to determine if save a sentence if it has no triples
def filter_out(data, flag=False):
    single,saved_data = [],[]
    for i, a_data in enumerate(data):
        sent_text = a_data['sentText']
        sent_words = nltk.word_tokenize(sent_text)
        triples_ = a_data['relationMentions']
        triples = set()
        for triple in triples_:
            if triple['label'] != 'N/A':
                triples.add((triple['em1Text'], triple['em2Text'], triple['label']))
        if len(sent_words) <= 200:
            if len(triples)==1:
                single.append(a_data)
            saved_data.append(a_data)
        elif flag:
            single.append(a_data)
        if (i + 1) * 1.0 % 10000 == 0:
            print('finish %f, %d/%d' % ((i + 1.0) / len(data), (i + 1), len(data)))
    print('instance number %d/%d' % (len(saved_data), len(data)))
    print('single triples number %d'%len(single))
    return saved_data

def write_data(f, data):
    out_data = [json.dumps(d,ensure_ascii=False) for d in data]
    for d in out_data:

        f.write(d)
        f.write('\n')
    f.close()

def run_filter():
    data = txt2json(test_file)
    saved_data = filter_out(data)
    f = open(train_filtered_filename, 'w',encoding='utf-8')
    write_data(f, saved_data)
    print('filter finish, saved in %s' % train_filtered_filename)

# run_filter()

def split(data):
    test_instance_num = 2000
    idx = random.sample(range(len(data)), test_instance_num)
    assert len(idx) == test_instance_num
    idx = set(idx)
    test_data = []
    train_data = []
    for i, a_data in enumerate(data):
        if i in idx:
            test_data.append(a_data)
        else:
            train_data.append(a_data)
    print(train_data[0])
    valid_instance_num = 2000
    valid_data = train_data[:valid_instance_num]
    train_data = train_data[valid_instance_num:]
    assert len(valid_data) == valid_instance_num
    assert len(test_data) == test_instance_num
    assert len(test_data) + len(train_data) + len(valid_data) == len(data)
    return test_data, train_data, valid_data

def run_split():
    data = read_json(train_filtered_filename)
    print('splitting')
    test_data, train_data, valid_data = split(data)
    print('saving')
    write_data(open(raw_test_filename, 'w',encoding='utf-8'), test_data)
    write_data(open(raw_train_filename, 'w',encoding='utf-8'), train_data)
    write_data(open(raw_valid_filename, 'w',encoding='utf-8'), valid_data)

# run_split()

def static_words(data):
    words = set()
    for a_data in data:
        sent_text = a_data['sentText']
        sent_words = nltk.word_tokenize(sent_text)
        words.update(set(sent_words))
        for rel in a_data['relationMentions']:
            em2Text=rel["em2Text"]
            em1Text = rel["em1Text"]
            words.add(em1Text)
            words.add(em2Text)
    words = list(words)
    words.insert(0, 'UNK')
    print('words number %d' % len(words))
    return list(words)

def static_relations(data):
    relations = set()
    for a_data in data:
        triples = a_data['relationMentions']
        for triple in triples:
            relation = triple['label']
            relations.add(relation)
    if 'N/A' in relations:
        relations.remove('N/A')
    relations = list(relations)
    # relations.insert(0, 'N/A')
    json.dump(relations, open(relations2id_filename, 'w'))
    print('relation number %d' % len(relations))
    return list(relations)

def run_static():
    data = read_json(train_filtered_filename)
    #decode_data=read_json(train_filtered_filename2)
    words = static_words(data+decode_data)
    words2id = dict()
    for i, w in enumerate(words):
        w=w.strip()
        words2id[w] = i
    json.dump(words2id, open(words2id_filename, 'w',encoding='utf-8'), indent=True)

    relations = static_relations(data)
    relations2id = dict()
    for i, r in enumerate(relations):
        relations2id[r] = i
    json.dump(relations2id, open(relations2id_filename, 'w',encoding='utf-8'), indent=True)

# run_static()


###################################
###      biobert的训练语料生成    ###
###################################
def triple_type(filename):
    data = read_json(filename)
    #single_triple=0
    multi_triple = 0
    tps,triples=[],[]
    for triple in data:
        relationMentions=triple["relationMentions"]
        if len(relationMentions)==1:
            continue
        else:
            for i in range(len(relationMentions)):
                multi_triple+=1
                sentText = triple["sentText"]
                em1Text=relationMentions[i]["em1Text"]
                em2Text = relationMentions[i]["em2Text"]
                label = relationMentions[i]["label"]
                if em1Text in sentText and em2Text in sentText:
                    sentText=sentText.replace(em1Text,"@GENE$")
                    sentText = sentText.replace(em2Text, "@DISEASE$")
                    triples.append((sentText,label))
                    tps.append(label)
    print(multi_triple,Counter(tps))
    return  triples
# triple_type(train_filtered_filename)


def write_tsv(f,data):
    with open(f, 'w',encoding='utf-8',newline='') as fw:
        tsv_w = csv.writer(fw, delimiter='\t')
        if "test" in f:
            tsv_w.writerow(['index', 'sentence', 'label'])
            idx=0
            for triple in data:
                text=triple[0]
                label="Expression"
                tsv_w.writerow([idx, text, label])
                idx+=1
        else:
            for triple in data:
                text=triple[0]
                label=triple[1]
                tsv_w.writerow([text, label])
def split_tsv():
    data = triple_type(train_filtered_filename)
    # print('splitting')
    # test_data, train_data, valid_data = split(data)
    print('saving')
    write_tsv(test_f_multi, data)



###################################
###      biobert的结果合并   ###
###################################

def write_tsv(f,data):
    with open(f, 'w',encoding='utf-8',newline='') as fw:
        tsv_w = csv.writer(fw, delimiter='\t')
        if "test" in f:
            tsv_w.writerow(['index', 'sentence', 'label'])
            idx=0
            for triple in data:
                text=triple[0]
                label="Expression"
                tsv_w.writerow([idx, text, label])
                idx+=1



def read_vec_bin():
    all_w2vec = dict()
    f = open(vector_bin_filename, 'r')
    for line in f:
        segs = line.strip().split(' ')
        word = segs[0]
        vector = [float(x) for x in segs[1:]]
        all_w2vec[word] = vector
    print('size %d' % len(all_w2vec))
    return all_w2vec

def load_words():
    return json.load(open(words2id_filename, 'r'))

def forword_vectors(words2id, all_w2vec):
    dim = len(all_w2vec[','])
    w2vec = dict()
    for w, idx in words2id.items():
        w2vec[idx] = all_w2vec.get(w, list(np.random.uniform(0, 1, dim))) #word没有对应的向量，就随机初始化
    assert len(w2vec) == len(words2id)
    return w2vec

def run_word_vectors():
    print('reading nyt_vec.bin')
    all_w2vec =read_vec_bin()
    words2id = load_words()
    print('prepare w2vec')
    w2vec = word_vectors(words2id, all_w2vec)
    print('dumping')
    json.dump(w2vec, open(words_id2vector_filename, 'w'))

def check(f):
    with open(f, 'r') as fr:
        count=[]
        c=0
        data = json.load(fr)
        all_sent_length, all_sent_id, all_triples_id = data
        for i,triple in enumerate(all_triples_id):
            count.append(len(triple))
            c+=1
        print(Counter(count))
        print(c)

def check_vector(f):
    with open(f,'r') as fr:
        count = []
        c = 0
        data = json.load(fr)
        print(len(list(data.values())[0]))

# check_vector(word2id)
# check_vector(id2vector)

def pubtator_api(pmid):
    url="https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/export/pubtator?pmids="+pmid
    response=requests.get(url)
    text=response.text.split("\n")
    abstract=text[0][len(pmid)+3:]+" "+text[1][len(pmid)+3:]
    return abstract


def response2txt():
    pubtator_home="E:/text_mining/0Data_Collection/211/"
    with open(pubtator_home+"pubtator_result","r",encoding="utf-8") as f,open(pubtator_home+"pubtator_result2","w",encoding="utf-8") as fw:
        count=0
        pmids=[]
        for line in f.readlines():
            if len(line.strip().split("\t"))==2:
               pmids.append(line.strip().split('\t')[0])
        print("there are %d abstracts loaded!" %len(pmids))
        for pmid in pubtator_pmid:
            if str(pmid) not in pmids:
                try:
                    abstract=pubtator_api(str(pmid))
                    fw.write(str(pmid)+"\t"+abstract+'\n')
                    count+=1
                    print("finishing getting %d abstracts"% count)
                except:
                    print(pmid)
        fw.close()