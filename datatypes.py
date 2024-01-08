from typing import Any, Optional
from dataclasses import dataclass
import pickle
import redis
import redis.client
import pandas as pd



null: Any = None # Lets be as un-Pythonic as possible to keep the neovim type linter happy 

class EntityObj:
    acc: Any

    def pkl(self):
        return pickle.dumps(self)

    @classmethod
    def unpkl(cls, pickled):
        return pickle.loads(pickled)

    def update(self, field, thing):
        if pd.isnull(thing):
            return
        # assert(hasattr(self, field)), (self, field, thing)
        if hasattr(self, field):
            ptr = self
        elif hasattr(self.acc, field):
            ptr = self.acc
        else:
            raise AttributeError(f"Object {self} has no attribute {field}")

        pre = getattr(ptr, field) 
        if pre == null:
            setattr(ptr, field, thing)
        elif isinstance(pre, set):
            setattr(ptr, field, pre | {thing})
        else:
            setattr(ptr, field, set([pre, thing]))


@dataclass
class GeneAtts:
    uniprot_kb_gene: Optional[str] = null 
    ensembl_gene: Optional[str] = null
    geneid: Optional[str] = null
    hgnc_id: Optional[str] = null
    hgnc_symbol: Optional[str] = null

@dataclass
class IsoformAtts:
    ensembl_trans: Optional[str] = null
    refseq_mrna: Optional[str] = null

@dataclass
class ProteinAtts:
    ensembl_prot: Optional[str] = null
    uniprot_swissprot: Optional[str] = null
    uniprot_trembl: Optional[str] = null

@dataclass
class ProteoformAtts:
    uniprot_isoform: Optional[str] = null
    refseq_prot: Optional[str] = null


class Gene(EntityObj):
    def __init__(self, identifier):
        self.identifier = identifier
        self.acc = GeneAtts()
        self.transcripts = set() 

class Isoform(EntityObj):
    def __init__(self, identifier):
        self.identifier = identifier
        self.acc = IsoformAtts()

        self.parent_gene = null
        self.proteins = set() 

class Protein(EntityObj):
    def __init__(self, identifier):
        self.identifier = identifier
        self.acc = ProteinAtts()

        self.parent_gene = null
        self.parent_isoform = null
        self.proteoforms = set() 

class Proteoform(EntityObj):
    def __init__(self, identifier):
        self.identifier = identifier
        self.acc = ProteoformAtts()

        self.parent_gene = null
        self.parent_isoform = null
        self.parent_protein = null


class ObjectyRedis(redis.Redis):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def oget(self, key):
        obj = self.get(key)
        if obj is None:
            return None 
        return pickle.loads(obj) # type: ignore

    def oset(self, key, obj):
        self.set(key, pickle.dumps(obj))

    def register_identifiers(self, obj):
        for key, val in obj.acc.__dict__.items():
            if key != 'identifier' and 'parent' not in key and val is not null and isinstance(val, str):
                keystr = 'KEY:' + val
                if self.exists(keystr):
                    pre = self.oget(keystr)
                    if pre == obj.identifier or (isinstance(pre, set) and obj.identifier in pre):
                        continue
                    else:
                        # print(f"Overloaded identifier: {keystr} ({pre}, {obj.identifier})")
                        if isinstance(pre, set):
                            self.oset(keystr, pre | {obj.identifier})
                        else:
                            self.oset(keystr, set([pre, obj.identifier]))
                else:
                    self.oset(keystr, obj.identifier)
        
        self.oset(obj.identifier, obj)

