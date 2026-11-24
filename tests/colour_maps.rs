//! The colour map shared by KoflCraft, FoehnCraft and MurmelCraft
//! (maintainer's decision, 2026-11-24): docs/colour-maps.md holds the eight
//! control points of ParaView's "Rainbow Desaturated", the same table the
//! other two projects write out and test. A changed point or position
//! breaks the match between the three projects' pictures here.

const POINTS: [(&str, &str); 8] = [
    ("0.0", "#4747DB"),
    ("0.143", "#00005C"),
    ("0.285", "#00FFFF"),
    ("0.429", "#008000"),
    ("0.571", "#FFFF00"),
    ("0.714", "#FF6100"),
    ("0.857", "#6B0000"),
    ("1.0", "#E04D4D"),
];

#[test]
fn colour_map_page_holds_the_shared_table() {
    let page = std::fs::read_to_string(
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("docs/colour-maps.md"),
    )
    .expect("docs/colour-maps.md");
    let rows: Vec<(String, String)> = page
        .lines()
        .filter_map(|l| {
            let cells: Vec<&str> = l
                .trim()
                .trim_matches('|')
                .split('|')
                .map(str::trim)
                .collect();
            match cells.as_slice() {
                [at, colour] if colour.starts_with("`#") => {
                    Some((at.to_string(), colour.trim_matches('`').to_string()))
                }
                _ => None,
            }
        })
        .collect();
    let want: Vec<(String, String)> = POINTS
        .iter()
        .map(|(a, c)| (a.to_string(), c.to_string()))
        .collect();
    assert_eq!(rows, want);
}
