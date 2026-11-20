//! v4zip library bindings / workspace placeholder
pub fn version() -> &'static str {
    "3.0.0"
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_version() {
        assert_eq!(version(), "3.0.0");
    }
}
