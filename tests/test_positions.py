import dataclasses

import gemmi
import pytest

from antibody_utils.numbering import Chain, Position, Scheme, parse_chain_position
from antibody_utils.numbering.positions import _CaseInsensitiveStrEnum


@pytest.mark.parametrize("insertion", ["", " ", "  ", "\t"])
def test_whitespace_means_no_insertion(insertion):
    assert Position(52, insertion) == Position(52)
    assert Position(52, insertion).insertion == ""


def test_positions_are_immutable():
    position = Position(52)
    with pytest.raises(dataclasses.FrozenInstanceError):
        position.insertion = " "  # type: ignore[misc]
    assert dataclasses.replace(position, insertion=" ").insertion == ""


@pytest.mark.parametrize("insertion", ["a", "AB", "1", "-"])
def test_invalid_insertion(insertion):
    with pytest.raises(ValueError, match="Insertion code"):
        Position(52, insertion)


def test_ordering():
    positions = [Position(36), Position(35, "B"), Position(35), Position(35, "A")]
    assert sorted(positions) == [
        Position(35),
        Position(35, "A"),
        Position(35, "B"),
        Position(36),
    ]


def test_hashable():
    assert {Position(52, "A"), Position(52, "A")} == {Position(52, "A")}


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        ("100", Position(100)),
        ("100 ", Position(100)),
        (" 100", Position(100)),
        ("100A", Position(100, "A")),
        ("100a", Position(100, "A")),
        (" 100 A ", Position(100, "A")),
    ],
)
def test_parse(text, expected):
    assert Position.parse(text) == expected


@pytest.mark.parametrize("text", ["", "A100", "100AB", "-1", "H100"])
def test_parse_invalid(text):
    with pytest.raises(ValueError, match="Not a valid position"):
        Position.parse(text)


def test_round_trip_through_str():
    for position in (Position(1), Position(111, "C")):
        assert Position.parse(str(position)) == position


def test_from_tuple():
    assert Position.from_tuple((111, "A")) == Position(111, "A")
    assert Position.from_tuple((111, " ")) == Position(111)
    assert Position.from_tuple(("111", " ")) == Position(111)


@pytest.mark.parametrize("position", [Position(1), Position(111, "C")])
def test_round_trip_through_tuple(position):
    assert Position.from_tuple(position.to_tuple()) == position


def test_to_tuple_matches_gemmi():
    for position in (Position(52), Position(52, "A")):
        seqid = gemmi.SeqId(*position.to_tuple())
        assert (seqid.num, seqid.icode) == position.to_tuple()


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        ("H100A", (Chain.HEAVY, Position(100, "A"))),
        ("l27", (Chain.LIGHT, Position(27))),
    ],
)
def test_parse_chain_position(text, expected):
    assert parse_chain_position(text) == expected


@pytest.mark.parametrize("text", ["100A", "K27", "H"])
def test_parse_chain_position_invalid(text):
    with pytest.raises(ValueError, match="chain-qualified"):
        parse_chain_position(text)


def test_enums_are_case_insensitive():
    assert Chain("l") is Chain.LIGHT
    assert Scheme("IMGT") is Scheme.IMGT
    with pytest.raises(ValueError, match="not a valid Scheme"):
        Scheme("wolfguy")


def test_case_insensitive_lookup_of_mixed_case_values():
    class Example(_CaseInsensitiveStrEnum):
        AHO = "AHo"

    for value in ("AHo", "aho", "AHO", "aHO"):
        assert Example(value) is Example.AHO
    with pytest.raises(ValueError, match="'ah' is not a valid"):
        Example("ah")
