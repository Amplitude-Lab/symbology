export const api = new Proxy({}, { get: () => async () => ({}) })
