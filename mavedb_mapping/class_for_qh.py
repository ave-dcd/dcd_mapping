from gene.database import AbstractDatabase
class ConcreteClass(AbstractDatabase):
    def __init__(self):
        pass

    def list_tables(self):
        pass

    def drop_db(self):
        pass

    def initialize_db(self):
        pass

    def get_source_metadata(self):
        pass

    def get_record_by_id(self):
        pass

    def get_refs_by_type(self):
        pass

    def get_all_concept_ids(self):
        pass

    def add_source_metadata(self):
        pass

    def add_record(self):
        pass

    def add_merged_record(self):
        pass

    def update_merge_ref(self):
        pass

    def delete_normalized_concepts(self):
        pass

    def delete_source(self):
        pass

    def complete_write_transaction(self):
        pass

    def close_connection(self):
        pass

    def load_from_remote(self):
        pass

    def export_db(self):
        pass
