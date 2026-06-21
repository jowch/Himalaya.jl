import { describe, it, expect } from 'vitest'
import { eventCombo, SHORTCUTS, matchShortcut } from '../src/print/shell/shortcuts'

const ev = (init: Partial<KeyboardEventInit>) => new KeyboardEvent('keydown', init)
const keysOf = (id: string) => Object.values(SHORTCUTS).find(d => d.id === id)?.keys ?? []

describe('shortcuts registry (rev-2 re-axis)', () => {
  it('CapsLock-X normalizes to x (existing behaviour)', () => {
    // Sanity-check the pre-existing normalization is still intact.
    expect(eventCombo(ev({ key: 'X' }))).toBe('x')
  })

  it('? normalizes to a stable token regardless of Shift', () => {
    expect(eventCombo(ev({ key: '?', shiftKey: true }))).toBe('?')
  })

  it('arrows are the sample/detail axis; [ and ] are gone', () => {
    expect(matchShortcut(ev({ key: 'ArrowUp' }))).toBe('prevSample')   // string id, not object
    expect(matchShortcut(ev({ key: 'ArrowLeft' }))).toBe('prevDetail')
    const allKeys = Object.values(SHORTCUTS).flatMap(d => d.keys)
    expect(allKeys).not.toContain('[')
    expect(allKeys).not.toContain(']')
  })

  it('new verbs are bound', () => {
    expect(keysOf('openFocus')).toContain('Enter')
    expect(keysOf('restore')).toContain('Backspace')
    expect(matchShortcut(ev({ key: ' ' }))).toBe('toggleSelect')
  })

  it('every combo resolves to exactly one id (no key bound twice)', () => {
    const seen = new Map<string, number>()
    for (const d of Object.values(SHORTCUTS)) for (const k of d.keys) seen.set(k, (seen.get(k) ?? 0) + 1)
    for (const [k, n] of seen) expect(`${k}:${n}`).toBe(`${k}:1`)
  })
})
