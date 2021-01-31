//! A logging implementation that emits status every "n" records that it is tracking.
use std::clone::Clone;
use std::cmp::PartialEq;
use std::default::Default;
use std::fmt::Debug;
use std::time::{Duration, SystemTime, SystemTimeError};

use case::CaseExt;
use chrono::{DateTime, Local};
use log::*;

/// By default, log every `n` records.
pub const DEFAULT_LOG_EVERY: usize = 100_000;

/// An interface for any struct that tracks and logs progress of records.
pub trait RecordLogger {
    /// Returns the elapsed time since this progress logger was created.
    fn elapsed_time(&self) -> Result<Duration, SystemTimeError>;
    /// Emit a log statement to log the most recent record added to the logger.
    fn emit(&mut self) -> Result<(), SystemTimeError>;
    /// Observe a record and optionally log it.
    fn observe(&mut self) -> Result<(), SystemTimeError>;
}

/// A progress logger that emits status every "n" records that it is tracking.
#[derive(Clone, PartialEq, Debug)]
pub struct ProgressLogger<'a> {
    /// The verb to use when emitting log statements, for example "processed", "read", "written.
    pub verb: &'a str,
    /// The noun of the record we are tracking, for example "reads" or "variants".
    pub noun: &'a str,
    /// How many records to collect before emitting a log.
    pub every: usize,
    count: usize,
    last_emit: usize,
    start: SystemTime,
    since: SystemTime,
}

impl<'a> ProgressLogger<'a> {
    /// Build a new empty progress logger, set the start time to now.
    pub fn new(verb: &'a str, noun: &'a str, every: usize) -> ProgressLogger<'a> {
        let now = SystemTime::now();
        let dt: DateTime<Local> = now.into();
        info!(
            "Logging records started at: {}.",
            dt.format("%Y-%m-%d %H:%M:%S")
        );
        ProgressLogger {
            verb,
            noun,
            every,
            count: 0,
            last_emit: 0,
            start: now,
            since: now,
        }
    }

    /// Format a duration for logging.
    fn format_duration(duration: Duration) -> String {
        let seconds = duration.as_secs() % 60;
        let minutes = (duration.as_secs() / 60) % 60;
        let hours = (duration.as_secs() / 60) / 60;
        format!("{:0>2}:{:0>2}:{:0>2}", hours, minutes, seconds)
    }
}

impl Default for ProgressLogger<'_> {
    /// Build a default empty progress logger, set the start time to now.
    fn default() -> Self {
        let now = SystemTime::now();
        ProgressLogger {
            verb: "processed",
            noun: "records",
            every: DEFAULT_LOG_EVERY,
            count: 0,
            last_emit: 0,
            start: now,
            since: now,
        }
    }
}

impl RecordLogger for ProgressLogger<'_> {
    /// Returns the elapsed time since this progress logger was created.
    fn elapsed_time(&self) -> Result<Duration, SystemTimeError> {
        SystemTime::now().duration_since(self.start)
    }

    /// Emit a log statement to log the most recent record added to the logger.
    fn emit(&mut self) -> Result<(), SystemTimeError> {
        let now = SystemTime::now();
        let since = ProgressLogger::format_duration(now.duration_since(self.since)?);
        let elapsed = ProgressLogger::format_duration(self.elapsed_time()?);
        info!(
            "{} {} {}. Elapsed time: {}. Time for last {} {}: {}",
            self.verb.to_capitalized(),
            self.count,
            self.noun,
            elapsed,
            self.count - self.last_emit,
            self.noun,
            since
        );
        self.last_emit = self.count;
        self.since = now;
        Ok(())
    }

    /// Observe a record and optionally log it.
    fn observe(&mut self) -> Result<(), SystemTimeError> {
        self.count += 1;
        if self.count % self.every == 0 {
            return self.emit();
        };
        Ok(())
    }
}
