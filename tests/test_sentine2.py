import pickle  # noqa
import random
from pathlib import Path

from satellitetools.common.sentinel2 import Sentinel2DataCollection

with open(
    Path(__file__).parent / "test_data" / "s2_data_collection_gee.pkl", "rb"
) as f:
    test_collection_gee: Sentinel2DataCollection = pickle.load(f)  # noqa


class TestSentinel2DataCollection:

    def test_filter_s2_items_by_tile(self):
        test_collection_gee.filter_s2_items_by_tile()
        assert len(test_collection_gee.s2_items) != 1
        test_collection_gee.filter_s2_items_by_tile("not_a_valid_tile")
        assert len(test_collection_gee.s2_items) == 0

    def test_sort_s2_items(self):
        original_s2_items = test_collection_gee.s2_items
        shuffled_s2_items = original_s2_items.copy()
        random.shuffle(shuffled_s2_items)
        test_collection_gee.s2_items = shuffled_s2_items
        test_collection_gee.sort_s2_items()
        assert test_collection_gee.s2_items == original_s2_items
