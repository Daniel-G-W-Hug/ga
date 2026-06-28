#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include <QFont>
#include <QString>

// ---------------------------------------------------------------------------
// Font sizing (DPI-independent — keeps the macOS-tuned layout on Windows too)
// ---------------------------------------------------------------------------
//
// All ga_view text (scene labels/legends, axis labels, the diagram-legend box,
// the status bar) is sized with setPixelSize() via these helpers rather than
// setPointSize() / the QFont point-size constructor.
//
// A point is 1/72 inch, so Qt converts a point size to pixels by the screen's
// *logical* DPI: 72 on macOS (1 pt -> 1 px) but 96 on Windows (1 pt -> 1.333
// px). All scene geometry, axis ticks and the legend box are laid out in
// pixels, so point-sized text comes out ~33 % larger on Windows relative to the
// drawing — breaking a layout tuned on macOS (text overflows the legend box,
// labels collide). A *pixel* size renders identically on every platform.
//
// The numeric size passed here is the pixel height the original macOS point
// size produced (pt == px at 72 dpi), so the macOS layout is preserved verbatim
// on Windows. Size text through these helpers, not raw setPointSize() / the
// QFont(family, pointSize, ...) constructor, to keep that invariant.

// fixed-width (monospace) font at a DPI-independent pixel size. Uses "Courier
// New" rather than the system fixed font (QFontDatabase::FixedFont) because the
// latter resolves to a different typeface per platform (Menlo on macOS, Consolas
// on Windows), so glyph widths drift between OSes. Courier New ships on both with
// matching metrics; TypeWriter is the fallback hint if it is ever absent.
inline QFont mono_font(int px, bool bold = false)
{
    QFont f("Courier New");
    f.setStyleHint(QFont::TypeWriter);
    f.setPixelSize(px);
    f.setBold(bold);
    return f;
}

// named proportional font (e.g. "Arial") at a DPI-independent pixel size.
// Prefer a family present on BOTH macOS and Windows (Arial, not Helvetica): a
// macOS-only family is silently substituted on Windows, shifting glyph metrics.
inline QFont named_font(char const* family, int px, int weight = QFont::Normal,
                        bool italic = false)
{
    QFont f(QString::fromLatin1(family), -1 /*point size unset*/, weight, italic);
    f.setPixelSize(px);
    return f;
}
