//! Memory allocator and memory-mapping checks.
//!
//! maligno uses mimalloc on every target. musl's built-in allocator (used by the
//! static Linux binary) splits the address space into tens of thousands of
//! separate mappings under the external sort, until Linux's per-process limit
//! (`vm.max_map_count`) is reached and allocations fail with RAM still free.
//! mimalloc keeps memory in a few large mappings.

use std::alloc::{GlobalAlloc, Layout};
use std::io::Write;
use std::sync::atomic::{AtomicBool, Ordering};

use mimalloc::MiMalloc;

/// The global allocator: mimalloc, plus a clear message on stderr if an
/// allocation fails (Rust then aborts as usual).
pub(crate) struct ReportingAlloc;

static REPORTED: AtomicBool = AtomicBool::new(false);

const ALLOC_FAILED_MSG: &[u8] = b"maligno: a memory allocation failed. Either the machine or a \
limit (ulimit -v, cgroup) ran out of memory, or the process reached Linux's limit on memory \
mappings (vm.max_map_count). Try a smaller --sort-mem, or pre-sort the inputs and use \
--presorted.\n";

/// Print the failure message once. The guard also stops recursion if writing
/// to stderr itself needs memory.
fn report_if_null(p: *mut u8) -> *mut u8 {
    if p.is_null() && !REPORTED.swap(true, Ordering::Relaxed) {
        let _ = std::io::stderr().write_all(ALLOC_FAILED_MSG);
    }
    p
}

unsafe impl GlobalAlloc for ReportingAlloc {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        report_if_null(unsafe { MiMalloc.alloc(layout) })
    }

    unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
        report_if_null(unsafe { MiMalloc.alloc_zeroed(layout) })
    }

    unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, new_size: usize) -> *mut u8 {
        report_if_null(unsafe { MiMalloc.realloc(ptr, layout, new_size) })
    }

    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        unsafe { MiMalloc.dealloc(ptr, layout) }
    }
}

/// True when `n` mappings are at least half of the `limit`.
fn at_least_half(n: u64, limit: u64) -> bool {
    limit > 0 && 2 * n >= limit
}

/// On Linux, warn when the process uses at least half of the system's
/// memory-mapping limit. Does nothing where `/proc` can't be read (e.g. macOS).
pub(crate) fn warn_if_many_memory_regions() {
    let (Ok(maps), Ok(limit)) = (
        std::fs::read_to_string("/proc/self/maps"),
        std::fs::read_to_string("/proc/sys/vm/max_map_count"),
    ) else {
        return;
    };
    let Ok(limit) = limit.trim().parse::<u64>() else {
        return;
    };
    let n = maps.lines().count() as u64;
    if at_least_half(n, limit) {
        eprintln!(
            "  WARNING: maligno is using {n} memory mappings, {}% of this system's limit \
             (vm.max_map_count = {limit}). If the limit is reached, memory allocations fail. \
             Use a smaller --sort-mem, or pre-sort the inputs and use --presorted.",
            n * 100 / limit
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn half_of_the_limit_triggers_the_warning() {
        assert!(!at_least_half(32_764, 65_530));
        assert!(at_least_half(32_765, 65_530));
        assert!(at_least_half(65_530, 65_530));
        assert!(!at_least_half(10, 0));
    }

    #[test]
    fn the_check_runs_without_proc() {
        warn_if_many_memory_regions();
    }
}
